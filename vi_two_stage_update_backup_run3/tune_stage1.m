%% TUNE_STAGE1  EO tuning of the fast inner stage (FO-PID + virtual rotor)
%  Design vector x = [Kp Ki Kd lam mu Kvi Dvi]
%  Stage 2 is off except for its minimum handover gain kh_min, which the
%  final controller always provides (see README, finding 1).
clear; clc; close all;

P   = mg_params();
S0  = stage2_defaults();
opt = sim_options();
opt.dt       = 2e-3;        % coarser step for speed (verify at 1e-3 afterwards)
opt.decim    = 50;
opt.kh_fixed = S0.kh_min;

ctx.P   = P;
ctx.opt = opt;
ctx.Tf  = 0.01;             % derivative filter; 1e-3 needs dt = 1e-3 (slower)
ctx.Tm_robust = 0.02;       % must stay stable with a 20 ms measurement lag
ctx.max_re    = -1e-4;      % all modes strictly stable
ctx.zeta_min  = 0.10;       % and every oscillatory mode damped >= 10 %
ctx.Tm_robust2 = 0.05;      % run 2: also stay well damped with a 50 ms lag
ctx.zeta_min2  = 0.05;      %        (run 1 gave zeta = 0.035 there)
% run 3: both gates must hold for every gain scaling stage 2 can apply
% (run 2 went unstable for alpha <= 0.1, which forced stage 2 to alpha = 1)
ctx.alpha_set  = [0.05 0.1 0.2 0.35 0.5 0.75 1 1.25 1.5];
ctx.w_rocof   = 0.2;
ctx.w_u       = 0.05;
ctx.scen = {mg_scenarios('I'), mg_scenarios('IV')};

% ---- normalisation: base-paper EO-FOPI with the same handover gain ----
Cref = ctrl_presets('eo_fopi');  Cref.Tf = ctx.Tf;
for s = 1:numel(ctx.scen)
    o = sim_mg(P, Cref, [], ctx.scen{s}, opt);
    ctx.ref(s) = struct('ITAE', max(o.J.ITAE,1e-9), ...
                        'ITAE_rocof', max(o.J.ITAE_rocof,1e-9), ...
                        'ISU', max(o.J.ISU,1e-9));
end

% ---- search space ----
%      Kp   Ki   Kd   lam  mu    Kvi  Dvi
lb = [ 0    0    0    0.1  0.05  0    0  ];
ub = [ 1000 1000 50   1.0  0.95  20   20 ];   % run 2: widened (run 1 hit Kp, Ki, Dvi limits)

eo.N = 40; eo.MaxIt = 60; eo.seed = 7;
eo.useParallel = ~isempty(ver('parallel'));   % uses parfor if the toolbox exists
tic
[xbest, Jbest, conv] = eo_optimize(@(x) cost_stage1(x, ctx), lb, ub, eo);
fprintf('\nStage-1 tuning finished in %.1f min, J = %.4f\n', toc/60, Jbest);

C1 = vec2ctrl1(xbest, ctx.Tf);
disp(C1)
[z, mr] = loop_margin(P, C1, 1, S0.kh_min, ctx.Tm_robust, opt);
fprintf('With Tm = %.0f ms: max Re = %.4g, min damping ratio = %.3f\n', ...
        1e3*ctx.Tm_robust, mr, z);
[z2, mr2] = loop_margin(P, C1, 1, S0.kh_min, ctx.Tm_robust2, opt);
fprintf('With Tm = %.0f ms: max Re = %.4g, min damping ratio = %.3f\n', ...
        1e3*ctx.Tm_robust2, mr2, z2);
zA = inf;  mA = -inf;
for a = ctx.alpha_set
    for Tm = [0 ctx.Tm_robust]
        [za, ma] = loop_margin(P, C1, a, S0.kh_min, Tm, opt);
        zA = min(zA, za);  mA = max(mA, ma);
    end
end
fprintf('Over alpha in [%.2f, %.2f], Tm = 0 and %.0f ms: max Re = %.4g, min damping ratio = %.3f\n', ...
        min(ctx.alpha_set), max(ctx.alpha_set), 1e3*ctx.Tm_robust, mA, zA);
save('stage1_result.mat', 'C1', 'xbest', 'Jbest', 'conv', 'ctx', 'lb', 'ub');

% ---- quick look ----
figure; semilogy(conv, 'LineWidth', 1.5); grid on
xlabel('Iteration'); ylabel('Best J'); title('Stage-1 EO convergence');

optF = sim_options(); optF.kh_fixed = S0.kh_min;
sc   = mg_scenarios('IV');
oRef = sim_mg(P, Cref, [], sc, optF);
o1   = sim_mg(P, C1,   [], sc, optF);
figure; plot(oRef.t, oRef.df, '--', o1.t, o1.df, 'LineWidth', 1.2); grid on
legend('EO-FOPI (paper gains) + handover', 'Tuned stage 1 + handover');
xlabel('Time (s)'); ylabel('\Deltaf (Hz)'); title('Scenario IV');
