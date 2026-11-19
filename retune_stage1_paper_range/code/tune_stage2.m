%% TUNE_STAGE2  EO tuning of the slow SOC-aware outer stage, stage 1 frozen
%  Design vector x = [s1 gap a_min a_em a_hi Kr kh_min dkh]
clear; clc; close all;
load('stage1_result.mat', 'C1');

P   = mg_params();
opt = sim_options();
opt.dt    = 5e-4;          % exact VI form needs <= 0.5 ms (see tune_stage1)
opt.decim = 400;

ctx.P   = P;
ctx.opt = opt;
ctx.C1  = C1;
ctx.S0  = stage2_defaults();
ctx.Tm_robust = 0.02;
ctx.max_re    = -1e-4;
ctx.zeta_min  = 0.10;
ctx.Tm_robust2 = 0.05;     % run 2: same second gate as stage 1
ctx.zeta_min2  = 0.05;
ctx.w_thr = 0.2;           % ESS energy throughput (battery wear proxy)
ctx.w_soc = 1.0;           % SOC-violation integral, p.u.*s outside [0.1, 0.9]
ctx.w_bms = 5.0;           % penalty if the BMS limit is reached
% SOC-stress scenarios + one normal-SOC scenario (stage 2 must not hurt it)
ctx.scen = {mg_scenarios('S5_lowSOC'), mg_scenarios('S6_highSOC'), mg_scenarios('IV')};

% ---- normalisation: stage 1 + fixed handover, no scheduler ----
optRef = opt;  optRef.kh_fixed = ctx.S0.kh_min;
for s = 1:numel(ctx.scen)
    o = sim_mg(P, C1, [], ctx.scen{s}, optRef);
    ctx.ref(s) = struct('ITAE', max(o.J.ITAE,1e-9), ...
                        'throughput_pus', max(o.J.throughput_pus,1e-9));
end

%      s1    gap   a_min a_em  a_hi  Kr    kh_min dkh
lb = [ 0.06  0.05  0.0   0.2   1.0   0.0   0.01   0.0 ];
ub = [ 0.30  0.20  1.0   1.0   1.5   1.0   0.20   1.5 ];   % run 2: a_min, dkh widened

eo.N = 30; eo.MaxIt = 50; eo.seed = 11;
eo.useParallel = ~isempty(ver('parallel'));
tic
[xbest, Jbest, conv] = eo_optimize(@(x) cost_stage2(x, ctx), lb, ub, eo);
fprintf('\nStage-2 tuning finished in %.1f min, J = %.4f\n', toc/60, Jbest);

S2 = vec2stage2(xbest, ctx.S0);
disp(S2)
save('stage2_result.mat', 'S2', 'xbest', 'Jbest', 'conv', 'ctx', 'lb', 'ub');

% ---- quick look: low-SOC scenario ----
optF = opt;  optF.decim = 40;  optFr = optF;  optFr.kh_fixed = ctx.S0.kh_min;
sc = mg_scenarios('S5_lowSOC');
oA = sim_mg(P, C1, [], sc, optFr);
oB = sim_mg(P, C1, S2, sc, optF);
figure;
subplot(3,1,1); plot(oA.t, oA.df, '--', oB.t, oB.df); grid on; ylabel('\Deltaf');
legend('Stage 1 + fixed handover', 'Two-stage'); title('S5: low SOC');
subplot(3,1,2); plot(oA.t, oA.soc, '--', oB.t, oB.soc); grid on; ylabel('SOC');
subplot(3,1,3); plot(oB.t, oB.a_dis, oB.t, oB.a_ch, oB.t, oB.kh); grid on
legend('\alpha_{dis}', '\alpha_{ch}', 'k_h'); xlabel('Time (s)');
