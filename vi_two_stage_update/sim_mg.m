function out = sim_mg(P, C, S, sc, opt)
%SIM_MG  Fixed-step simulation of the islanded MG with the two-stage VI controller.
%
%   out = sim_mg(P, C, S, sc, opt)
%     P   : plant parameters               (mg_params)
%     C   : stage-1 controller             (ctrl_presets / vec2ctrl1)
%     S   : stage-2 scheduler parameters   (stage2_defaults) or [] = stage 2 off
%     sc  : scenario                       (mg_scenarios)
%     opt : numerical options              (sim_options)
%
%   Every linear block is discretised exactly (ZOH); GRC, valve limits, ESS
%   limits and the BMS SOC taper are applied between blocks each step.
%   Sign convention: dPVI > 0 means the ESS DISCHARGES into the grid.
%
%   out.t, out.df, out.Pvi, out.soc, out.Pm, out.a_dis, out.a_ch, out.kh, out.prec
%   out.J        : metrics computed at full resolution
%   out.unstable : true if the run diverged

if nargin < 5 || isempty(opt), opt = sim_options(); end
assert(C.lam > 0 && C.lam <= 1, 'sim_mg: lam must be in (0,1]');
assert(C.mu >= 0 && C.mu < 1,  'sim_mg: mu must be in [0,1)');

dt = opt.dt;
n  = round(sc.T / dt);
t  = (0:n-1)' * dt;
Dist = sc.dist(t);                         % n x 4 : [Presi Pind Pwind Psolar]
PL   = Dist(:,1) + Dist(:,2);
Pw0  = Dist(:,3);
Ppv0 = Dist(:,4);

% ---------- discretised first-order blocks ----------
ag  = exp(-dt/P.Tg);   at  = exp(-dt/P.Tt);
aw  = exp(-dt/P.Twt);  apv = exp(-dt/P.Tpv);
ae  = exp(-dt/P.Tinv);
asw = exp(-P.D*dt/(2*P.H));                % swing equation 1/(2Hs+D)
% Tf = 0 : EXACT virtual rotor + ESS, (Kvi*s + Dvi)/(1 + s*Tinv), realised
%          as in the base paper / Step-1 model (no numerical derivative).
% Tf > 0 : old form, Kvi*s/(Tf*s+1) followed by the ESS lag. NOTE: with
%          Tf = dt the ZOH filter over-states the derivative by
%          1/(1-exp(-1)) = 1.58x, so do not use Tf ~ dt.
Tf    = C.Tf;
exact = (Tf == 0);
if exact
    af = 0;
else
    if Tf < dt, warning('sim_mg: Tf < dt, derivative is poorly resolved'); end
    af = exp(-dt/Tf);
end
% OSTIT2F-FOPI (base paper): IT2 fuzzy self-tuning of Kp, Ki each step,
% 'operating point' form of the Step-1 model: Kp = Kp0 + (yP - UP/2)
isFz = isfield(C, 'type') && strcmp(C.type, 'it2f');
Kp = C.Kp;  Ki = C.Ki;  meas_prev = 0;
afa = exp(-dt/1.0);                        % 1 s filter on |df| for stage 2
if opt.Tm > 0, am = exp(-dt/opt.Tm); else, am = 0; end

% ---------- fractional operators ----------
% s^(-lam) = (1/s) * s^(1-lam) ; s^(mu) direct
[AI,BI,CI,DI] = oustaloup_ss(1 - C.lam, opt.wb, opt.wh, opt.N);
[AId,BId]     = zoh_c2d(AI, BI, dt);
[AD,BD,CD,DD] = oustaloup_ss(C.mu, opt.wb, opt.wh, opt.N);
[ADd,BDd]     = zoh_c2d(AD, BD, dt);
xI = zeros(size(AI,1),1);
xD = zeros(size(AD,1),1);
useI = (C.Ki ~= 0);
useD = (C.Kd ~= 0);

% ---------- measurement delay buffer ----------
nd  = round(opt.Td / dt);
buf = zeros(nd + 1, 1);
ptr = 1;

% ---------- ESS energy ----------
E = P.E_MWh;
if isfield(sc, 'E_MWh') && ~isempty(sc.E_MWh), E = sc.E_MWh; end
Epu = E * 3600 / P.Sbase;                  % p.u. * s

% ---------- states ----------
df = 0; xg = 0; yr = 0; xt = 0; xc = 0; xw = 0; xpv = 0; xe = 0;
xi = 0; xdf = 0; fa = 0; xm = 0;
wv = 0; xp = 0;                            % exact form: lag of u1 and of P_rec
soc = sc.soc0;

% ---------- stage 2 ----------
hasS  = ~isempty(S);
a_dis = 1; a_ch = 1; prec = 0; kh = opt.kh_fixed;
if hasS
    M2 = max(round(S.Ts2 / dt), 1);
    ka = exp(-S.Ts2 / S.tau_a);
    kr = exp(-S.Ts2 / S.tau_r);
    kh = S.kh_min;
end

% ---------- metric accumulators ----------
IAE = 0; ITAE = 0; ISE = 0; ITAEr = 0; maxRo = 0; OS = -inf; US = inf;
THR = 0; ISU = 0; socmin = soc; socmax = soc; socviol = 0; tout = 0;
unstable = false;

% ---------- storage ----------
dec = max(opt.decim, 1);
m   = floor((n - 1)/dec) + 1;
rec = zeros(m, 9);
j   = 0;

for k = 1:n
    tk = t(k);

    % ===== STAGE 2 (slow): SOC-aware scheduler =====
    if hasS && mod(k - 1, M2) == 0
        [ad_t, ac_t, kh_t, pr_t] = soc_scheduler(soc, fa, S);
        a_dis = ka*a_dis + (1 - ka)*ad_t;
        a_ch  = ka*a_ch  + (1 - ka)*ac_t;
        kh    = ka*kh    + (1 - ka)*kh_t;
        prec  = kr*prec  + (1 - kr)*pr_t;
    end

    % ===== measurement (lag + delay) =====
    if am > 0, meas = xm; else, meas = df; end
    if nd > 0
        buf(ptr) = meas;
        ptr = mod(ptr, nd + 1) + 1;
        meas = buf(ptr);
    end

    % ===== STAGE 1 (fast): FO-PID + virtual rotor =====
    e = -meas;
    if isFz
        [yP, yIt] = it2flc_gains(C.KE*meas, C.KDE*(meas - meas_prev)/dt, C.UP, C.UI);
        Kp = C.Kp0 + (yP  - C.UP/2);
        Ki = C.Ki0 + (yIt - C.UI/2);
        meas_prev = meas;
    end
    if useI, yI = CI*xI + DI*xi; else, yI = 0; end
    if useD, yD = CD*xD + DD*e;  else, yD = 0; end
    u1  = Kp*e + Ki*yI + C.Kd*yD;
    if exact
        Pc = C.Kvi*(u1 - wv)/P.Tinv + C.Dvi*wv;      % already includes the ESS lag
        if Pc >= 0, Pc = a_dis*Pc; else, Pc = a_ch*Pc; end
        xe = Pc + xp;
    else
        der = (u1 - xdf) / Tf;
        Pc  = C.Kvi*der + C.Dvi*u1;
        if Pc >= 0, Pc = a_dis*Pc; else, Pc = a_ch*Pc; end
        Pcmd = Pc + prec;
    end

    % ===== ESS output with power limits and BMS SOC taper =====
    c_dis = min(max((soc - P.soc_hmin)/P.soc_ramp, 0), 1);
    c_ch  = min(max((P.soc_hmax - soc)/P.soc_ramp, 0), 1);
    Pvi   = min(max(xe, P.Pmin*c_ch), P.Pmax*c_dis);

    % ===== plant =====
    Pnet = xt + xw + xpv + Pvi - PL(k);
    dfn  = asw*df + (1 - asw)*Pnet/P.D;
    xgn  = ag*xg + (1 - ag)*(xc - df/P.R);
    yrn  = yr + min(max(xg - yr, -P.GRC*dt), P.GRC*dt);     % GRC
    v    = min(max(yr, P.VL), P.VU);                         % valve limits
    xtn  = at*xt + (1 - at)*v;
    % secondary control (+ stage-2 handover of the VI share of ESS power)
    xc   = xc + dt*(-P.Ks*P.beta*df + kh*(Pvi - prec));
    xw   = aw*xw   + (1 - aw)*Pw0(k);
    xpv  = apv*xpv + (1 - apv)*Ppv0(k);
    if exact
        wv = ae*wv + (1 - ae)*u1;
        xp = ae*xp + (1 - ae)*prec;
    else
        xe = ae*xe + (1 - ae)*Pcmd;
    end
    if useI, xI = AId*xI + BId*xi; end
    xi   = xi + dt*e;
    if useD, xD = ADd*xD + BDd*e; end
    xdf  = af*xdf + (1 - af)*u1;
    fa   = afa*fa + (1 - afa)*abs(meas);
    if am > 0, xm = am*xm + (1 - am)*df; end
    if Pvi > 0, pe = Pvi/P.eta; else, pe = Pvi*P.eta; end
    soc  = soc - dt*pe/Epu;

    % ===== metrics =====
    rocof = (dfn - df)/dt;
    ad    = abs(df);
    IAE   = IAE  + ad*dt;
    ITAE  = ITAE + tk*ad*dt;
    ISE   = ISE  + df*df*dt;
    ITAEr = ITAEr + tk*abs(rocof)*dt;
    if abs(rocof) > maxRo, maxRo = abs(rocof); end
    if df > OS, OS = df; end
    if df < US, US = df; end
    THR   = THR + abs(Pvi)*dt;
    ISU   = ISU + Pvi*Pvi*dt;
    if soc < socmin, socmin = soc; end
    if soc > socmax, socmax = soc; end
    socviol = socviol + (max(P.soc_lo_ok - soc, 0) + max(soc - P.soc_hi_ok, 0))*dt;
    if soc < 0.2 || soc > 0.8, tout = tout + dt; end

    % ===== store =====
    if mod(k - 1, dec) == 0
        j = j + 1;
        rec(j,:) = [tk, df, Pvi, soc, xt, a_dis, a_ch, kh, prec];
    end

    % ===== advance =====
    df = dfn; xg = xgn; yr = yrn; xt = xtn;
    if ~isfinite(df) || abs(df) > 10
        unstable = true;
        break
    end
end

rec = rec(1:j,:);
out.t     = rec(:,1);  out.df    = rec(:,2);  out.Pvi = rec(:,3);
out.soc   = rec(:,4);  out.Pm    = rec(:,5);
out.a_dis = rec(:,6);  out.a_ch  = rec(:,7);  out.kh  = rec(:,8);
out.prec  = rec(:,9);
out.unstable = unstable;
out.J = struct('OS',OS,'US',US,'IAE',IAE,'ITAE',ITAE,'ISE',ISE, ...
               'ITAE_rocof',ITAEr,'maxRoCoF',maxRo, ...
               'throughput_pus',THR,'throughput_MWh',THR*P.Sbase/3600, ...
               'ISU',ISU,'socmin',socmin,'socmax',socmax,'socend',soc, ...
               'soc_viol',socviol,'t_outside_20_80',tout);
end
