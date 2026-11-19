function [ev, A] = mg_linearize(P, C, alpha, kh, Tm, opt)
%MG_LINEARIZE  Small-signal model with stage-2 outputs frozen.
%   Ignores GRC, valve/ESS limits and the BMS taper (valid around an
%   operating point inside all limits). alpha = frozen gain scaling,
%   kh = frozen handover gain, Tm = measurement lag (0 = none).
%   The fractional operators use the same Oustaloup approximation as
%   sim_mg, so integer-order eigenvalues are the right stability test.
%   States: [df xg xt xc xe xi xdf (xm) xI xD]
if nargin < 6 || isempty(opt), opt = sim_options(); end
[AI,BI,CI,DI] = oustaloup_ss(1 - C.lam, opt.wb, opt.wh, opt.N);
[AD,BD,CD,DD] = oustaloup_ss(C.mu,      opt.wb, opt.wh, opt.N);
nI   = size(AI,1);
nD   = size(AD,1);
hasM = Tm > 0;
n0   = 7 + hasM;
n    = n0 + nI + nD;
iI   = n0 + (1:nI);
iD   = n0 + nI + (1:nD);

A = zeros(n);
for jj = 1:n
    x = zeros(n,1); x(jj) = 1;
    A(:,jj) = f(x);
end
% drop states that cannot affect the loop
keep = true(n,1);
if C.Ki == 0, keep([6, iI]) = false; end
if C.Kd == 0, keep(iD) = false; end
if C.Tf == 0, keep(7) = false; end           % exact form has no derivative filter
A  = A(keep, keep);
ev = eig(A);

    function d = f(x)
        df = x(1); xg = x(2); xt = x(3); xc = x(4); xe = x(5);
        xi = x(6); xdf = x(7);
        if hasM, m = x(8); else, m = df; end
        xI = x(iI); xD = x(iD);
        e  = -m;
        yI = CI*xI + DI*xi;
        yD = CD*xD + DD*e;
        u1 = C.Kp*e + C.Ki*yI + C.Kd*yD;
        d = zeros(n,1);
        if C.Tf == 0
            % exact (Kvi*s + Dvi)/(1 + s*Tinv): state 5 = lag of u1 (wv),
            % ESS power is algebraic, state 7 is unused
            wv = xe;
            xe = alpha*(C.Kvi*(u1 - wv)/P.Tinv + C.Dvi*wv);
            Pc = u1;  der = 0;               % so that d(5) = (u1 - wv)/Tinv
        else
            der = (u1 - xdf)/C.Tf;
            Pc  = alpha*(C.Kvi*der + C.Dvi*u1);
        end
        d(1) = (xt + xe - P.D*df)/(2*P.H);
        d(2) = (xc - df/P.R - xg)/P.Tg;
        d(3) = (xg - xt)/P.Tt;
        d(4) = -P.Ks*P.beta*df + kh*xe;
        if C.Tf == 0, d(5) = (Pc - wv)/P.Tinv; else, d(5) = (Pc - xe)/P.Tinv; end
        d(6) = e;
        d(7) = der;
        if hasM, d(8) = (df - m)/Tm; end
        d(iI) = AI*xI + BI*xi;
        d(iD) = AD*xD + BD*e;
    end
end
