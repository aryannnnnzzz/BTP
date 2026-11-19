%% RUN_VALIDATION  Step 6: full comparison and robustness suite (dt = 0.5 ms).
%  Controllers (ablation order):
%   A  EO-FOPI, paper gains, no handover           (base-paper benchmark)
%   O  OSTIT2F-FOPI, paper Table 5, no handover    (base-paper proposed)
%   B  EO-FOPI + fixed handover kh_min             (effect of handover alone)
%   C  Tuned stage 1 + fixed handover              (effect of the new inner stage)
%   D  Two-stage: stage 1 + SOC-aware stage 2      (proposed)
%  All controllers use the exact VI block (Kvi*s + Dvi)/(1 + s*Tinv) (Tf = 0),
%  as in the base paper and the Step-1 model.
clear; clc; close all;
P  = mg_params();
load('stage1_result.mat', 'C1');  load('stage2_result.mat', 'S2');
S0 = stage2_defaults();

base = sim_options();  base.dt = 5e-4;  base.decim = 40;
Cp   = ctrl_presets('eo_fopi');  Cp.Tf = 0;
Co   = ctrl_presets('ostit2f');  Co.Tf = 0;

ctrl = struct('name', {}, 'C', {}, 'S', {}, 'kh', {});
ctrl(1).name = 'A: EO-FOPI (paper)';       ctrl(1).C = Cp; ctrl(1).S = []; ctrl(1).kh = 0;
ctrl(2).name = 'O: OSTIT2F-FOPI (paper)';  ctrl(2).C = Co; ctrl(2).S = []; ctrl(2).kh = 0;
ctrl(3).name = 'B: EO-FOPI + handover';    ctrl(3).C = Cp; ctrl(3).S = []; ctrl(3).kh = S0.kh_min;
ctrl(4).name = 'C: Stage 1 + handover';    ctrl(4).C = C1; ctrl(4).S = []; ctrl(4).kh = S2.kh_min;
ctrl(5).name = 'D: Two-stage (proposed)';  ctrl(5).C = C1; ctrl(5).S = S2; ctrl(5).kh = 0;

scen = {'I','II','III','IV','S5_lowSOC','S6_highSOC','S7_long'};
rows = {};
for s = 1:numel(scen)
    sc = mg_scenarios(scen{s});
    for c = 1:numel(ctrl)
        o = base;  o.kh_fixed = ctrl(c).kh;
        r = sim_mg(P, ctrl(c).C, ctrl(c).S, sc, o);
        res(s,c) = r; %#ok<SAGROW>
        rows(end+1,:) = {scen{s}, ctrl(c).name, r.J.OS, r.J.US, r.J.IAE, r.J.ITAE, ...
            r.J.maxRoCoF, r.J.socmin, r.J.socmax, r.J.socend, ...
            r.J.throughput_MWh, r.J.t_outside_20_80, r.unstable}; %#ok<SAGROW>
    end
end
T = cell2table(rows, 'VariableNames', {'Scenario','Controller','OS','US','IAE','ITAE', ...
    'maxRoCoF','SOCmin','SOCmax','SOCend','Throughput_MWh','t_out_20_80_s','Unstable'});
disp(T);  writetable(T, 'validation_results.csv');

% ---------- figures: frequency (all controllers) + SOC (proposed only) ----------
% The base paper does not model SOC, so SOC is shown only for the proposed
% two-stage controller (its SOC numbers for every controller stay in the table).
% Controllers whose df traces coincide (e.g. A and O, or C and D at normal
% SOC, where stage 2 stays idle) are drawn ONCE with a combined legend entry,
% so every legend entry has a visible line. Same colour per controller in
% every figure.
iD  = find(strncmp({ctrl.name}, 'D:', 2), 1);
co  = lines(numel(ctrl));
tol = 0.01;                                   % "identical" = within 1 % of the peak |df|
for s = 1:numel(scen)
    rep = 1:numel(ctrl);                      % rep(c) = controller whose line represents c
    for c = 2:numel(ctrl)
        for k = 1:c-1
            if rep(k) == k && numel(res(s,c).df) == numel(res(s,k).df) && ...
               max(abs(res(s,c).df - res(s,k).df)) <= tol*max(abs(res(s,k).df))
                rep(c) = k;  break
            end
        end
    end
    shown = find(rep == 1:numel(ctrl));
    pk    = arrayfun(@(c) max(abs(res(s,c).df)), shown);
    % a curve whose peak is < 10 % of the largest one is flat at full scale
    % (e.g. B and D next to the battery-limit spikes in S5/S6): show those in
    % a separate zoomed panel, and list each curve only where it is visible
    big   = pk >= 0.1*max(pk);
    if all(big), groups = {shown};  else, groups = {shown(big), shown(~big)}; end
    nr = numel(groups) + 1;
    figure('Name', scen{s});
    if nr == 3, set(gcf, 'Position', get(gcf, 'Position').*[1 1 1 1.35]); end
    for g = 1:numel(groups)
        subplot(nr,1,g); hold on; grid on
        h = [];  lab = {};
        % lines drawn first are thicker, so a nearly coinciding curve drawn
        % later (e.g. B over A at normal SOC) leaves the earlier one visible
        % as a coloured edge
        lw = max(2.6 - 0.8*(0:numel(groups{g})-1), 1.0);
        for ic = 1:numel(groups{g})
            c = groups{g}(ic);
            h(end+1) = plot(res(s,c).t, res(s,c).df, 'LineWidth', lw(ic), 'Color', co(c,:)); %#ok<AGROW>
            lab{end+1} = strjoin({ctrl(rep == c).name}, '  =  '); %#ok<AGROW>
        end
        ylabel('\Deltaf (Hz)');
        if g == 1
            title(['Scenario ' scen{s}], 'Interpreter', 'none');
        else
            title('Same scenario, zoomed on the small deviations', 'Interpreter', 'none');
        end
        legend(h, lab, 'Location', 'best', 'Interpreter', 'none');
    end
    cD = co(rep(iD), :);                      % colour of the line that shows D
    subplot(nr,1,nr); grid on
    plot(res(s,iD).t, res(s,iD).soc, 'LineWidth', 1.1, 'Color', cD);
    ylabel('SOC'); xlabel('Time (s)');
    title('SOC of the proposed two-stage controller');
end

% ---------- robustness: +/-20 % parameters (Scenario IV) ----------
fields = {'H','D','Tg','Tt','Twt','Tpv','R'};
fac    = [0.8 1.0 1.2];
sc4    = mg_scenarios('IV');
fprintf('\nParameter sensitivity, Scenario IV, two-stage (IAE):\n');
for f = 1:numel(fields)
    v = zeros(1, numel(fac));
    for k = 1:numel(fac)
        Pp = P;  Pp.(fields{f}) = P.(fields{f})*fac(k);
        r  = sim_mg(Pp, C1, S2, sc4, base);
        v(k) = r.J.IAE;
    end
    fprintf('  %-4s  -20%%: %.5f   nominal: %.5f   +20%%: %.5f\n', fields{f}, v);
end
% low-inertia stress
Pp = P;  Pp.H = 0.5*P.H;
r  = sim_mg(Pp, C1, S2, sc4, base);
fprintf('  H x 0.5: IAE %.5f, US %.5f\n', r.J.IAE, r.J.US);

% ---------- measurement lag / delay (Scenario IV) ----------
% dt = 0.2 ms here: with a lag, the ~20 Hz inertia mode is lightly damped and
% the 1 ms block-by-block discretisation adds enough phase lag to fake an
% oscillation. 0.2 ms results agree with the eigenvalue analysis.
fine = base;  fine.dt = 2e-4;  fine.decim = 100;
% sim_mg's "unstable" flag only catches divergence (|df| > 10 Hz). A loop with
% unstable eigenvalues is usually held in a bounded limit cycle by the ESS
% power limit, so also report: LC = 1 if the RMS of df over the last 10 s is
% above 1e-3 Hz (Scenario IV is quiet by then), and the linear max Re(eig).
lc = @(r) double(r.unstable || sqrt(mean(r.df(r.t > r.t(end) - 10).^2)) > 1e-3);
fprintf('\nMeasurement effects, Scenario IV, dt = 0.2 ms (IAE, diverged, limit cycle, max Re(eig), zeta):\n');
for Tm = [0 0.01 0.02 0.05]
    o = fine;  o.Tm = Tm;
    rp = sim_mg(P, Cp, [], sc4, o);
    ro = sim_mg(P, Co, [], sc4, o);
    rd = sim_mg(P, C1, S2, sc4, o);
    [zp, mrp] = loop_margin(P, Cp, 1, 0,         Tm, base);   % paper: no handover
    [zd, mrd] = loop_margin(P, C1, 1, S2.kh_min, Tm, base);
    fprintf(['  Tm = %2.0f ms | EO-FOPI paper: %9.5f (div %d, LC %d, Re %+8.3g, zeta %6.3f)' ...
             ' | OSTIT2F paper: %9.5f (div %d, LC %d)' ...
             ' | two-stage: %9.5f (div %d, LC %d, Re %+8.3g, zeta %6.3f)\n'], ...
            1e3*Tm, rp.J.IAE, rp.unstable, lc(rp), mrp, zp, ...
                    ro.J.IAE, ro.unstable, lc(ro), ...
                    rd.J.IAE, rd.unstable, lc(rd), mrd, zd);
end

% ---------- ESS time constant ----------
% With the exact VI form the direct df -> ESS gain is G = Kvi*(Kp+Kd*wh^mu)/Tinv,
% so a small Tinv needs a much smaller step: dt is chosen so dt*G/(2H) <= 0.2.
fprintf('\nESS time constant, Scenario IV (IAE, limit cycle, max Re(eig), zeta):\n');
for Tinv = [10 1 0.2]
    Pp = P;  Pp.Tinv = Tinv;
    G  = max(S2.a_hi, 1)*C1.Kvi*(C1.Kp + C1.Kd*base.wh^C1.mu)/Tinv;
    o  = fine;  o.dt = min(fine.dt, 0.2*2*P.H/G);  o.decim = max(1, round(0.01/o.dt));
    rd = sim_mg(Pp, C1, S2, sc4, o);
    [zd, mrd] = loop_margin(Pp, C1, 1, S2.kh_min, 0, base);
    fprintf('  Tinv = %4.1f s (dt = %.1e s) | two-stage IAE %.5f, US %.5f, LC %d, Re %+8.3g, zeta %.3f\n', ...
            Tinv, o.dt, rd.J.IAE, rd.J.US, lc(rd), mrd, zd);
end
