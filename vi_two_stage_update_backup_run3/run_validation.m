%% RUN_VALIDATION  Step 6: full comparison and robustness suite (dt = 1 ms).
%  Controllers (ablation order):
%   A  EO-FOPI, paper gains, no handover           (base-paper behaviour)
%   B  EO-FOPI + fixed handover kh_min             (effect of handover alone)
%   C  Tuned stage 1 + fixed handover              (effect of the new inner stage)
%   D  Two-stage: stage 1 + SOC-aware stage 2      (proposed)
%  OSTIT2F-FOPI is not re-implemented here: take its numbers from your
%  Step-1 Simulink model for Scenarios I-IV.
clear; clc; close all;
P  = mg_params();
load('stage1_result.mat', 'C1');  load('stage2_result.mat', 'S2');
S0 = stage2_defaults();

base = sim_options();  base.dt = 1e-3;  base.decim = 20;
Cp   = ctrl_presets('eo_fopi');  Cp.Tf = C1.Tf;

ctrl(1) = struct('name','A: EO-FOPI (paper)',     'C',Cp, 'S',[], 'kh',0);
ctrl(2) = struct('name','B: EO-FOPI + handover',  'C',Cp, 'S',[], 'kh',S0.kh_min);
ctrl(3) = struct('name','C: Stage 1 + handover',  'C',C1, 'S',[], 'kh',S2.kh_min);
ctrl(4) = struct('name','D: Two-stage (proposed)','C',C1, 'S',S2, 'kh',0);

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

% ---------- figures: frequency + SOC per scenario ----------
for s = 1:numel(scen)
    figure('Name', scen{s});
    subplot(2,1,1); hold on; grid on
    for c = 1:numel(ctrl), plot(res(s,c).t, res(s,c).df, 'LineWidth', 1.1); end
    ylabel('\Deltaf (Hz)'); title(['Scenario ' scen{s}], 'Interpreter', 'none');
    legend({ctrl.name}, 'Location', 'best');
    subplot(2,1,2); hold on; grid on
    for c = 1:numel(ctrl), plot(res(s,c).t, res(s,c).soc, 'LineWidth', 1.1); end
    ylabel('SOC'); xlabel('Time (s)');
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
    rd = sim_mg(P, C1, S2, sc4, o);
    [zp, mrp] = loop_margin(P, Cp, 1, S0.kh_min, Tm, base);
    [zd, mrd] = loop_margin(P, C1, 1, S2.kh_min, Tm, base);
    fprintf(['  Tm = %2.0f ms | EO-FOPI paper: %9.5f (div %d, LC %d, Re %+8.3g, zeta %6.3f)' ...
             ' | two-stage: %9.5f (div %d, LC %d, Re %+8.3g, zeta %6.3f)\n'], ...
            1e3*Tm, rp.J.IAE, rp.unstable, lc(rp), mrp, zp, ...
                    rd.J.IAE, rd.unstable, lc(rd), mrd, zd);
end

% ---------- ESS time constant ----------
fprintf('\nESS time constant, Scenario IV, dt = 0.2 ms (IAE, limit cycle, max Re(eig), zeta):\n');
for Tinv = [10 1 0.2]
    Pp = P;  Pp.Tinv = Tinv;
    rd = sim_mg(Pp, C1, S2, sc4, fine);
    [zd, mrd] = loop_margin(Pp, C1, 1, S2.kh_min, 0, base);
    fprintf('  Tinv = %4.1f s | two-stage IAE %.5f, US %.5f, LC %d, Re %+8.3g, zeta %.3f\n', ...
            Tinv, rd.J.IAE, rd.J.US, lc(rd), mrd, zd);
end
