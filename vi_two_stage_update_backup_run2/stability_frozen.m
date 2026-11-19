%% STABILITY_FROZEN  Step 8: frozen-gain stability map of the two-stage loop.
%  Because stage 2 is slow (tau_a = 2 s, Ts2 = 0.1 s) its outputs are treated
%  as constant parameters; the closed loop must be stable for every frozen
%  (alpha, kh) the scheduler can produce, with and without measurement lag.
clear; clc; close all;
P  = mg_params();  opt = sim_options();
load('stage1_result.mat', 'C1');  load('stage2_result.mat', 'S2');
Cpaper = ctrl_presets('eo_fopi');  Cpaper.Tf = C1.Tf;

alphas = linspace(min(S2.a_min, 0.05), S2.a_hi, 25);
khs    = linspace(S2.kh_min, S2.kh_max, 15);
Tms    = [0 0.01 0.02 0.05];

fprintf('Worst case over the (alpha, kh) grid: max Re(eig) and min damping ratio\n');
fprintf('%8s | %10s %8s | %10s %8s\n', 'Tm (ms)', 'two-stage', 'zeta', 'EO paper', 'zeta');
MR = zeros(numel(alphas), numel(khs));
for Tm = Tms
    worst = -inf;  worstP = -inf;  zw = inf;  zwP = inf;
    for ia = 1:numel(alphas)
        for ik = 1:numel(khs)
            [z, mr] = loop_margin(P, C1, alphas(ia), khs(ik), Tm, opt);
            MR(ia,ik) = mr;  worst = max(worst, mr);  zw = min(zw, z);
            [zp, mrp] = loop_margin(P, Cpaper, alphas(ia), khs(ik), Tm, opt);
            worstP = max(worstP, mrp);  zwP = min(zwP, zp);
        end
    end
    fprintf('%8.0f | %10.4g %8.3f | %10.4g %8.3f\n', 1e3*Tm, worst, zw, worstP, zwP);
    if Tm == 0.02
        figure; imagesc(khs, alphas, MR); axis xy; colorbar
        xlabel('k_h (1/s)'); ylabel('\alpha'); title('max Re(\lambda), T_m = 20 ms');
    end
end

% eigenvalue map at the nominal point
ev = mg_linearize(P, C1, 1, S2.kh_min, 0.02, opt);
figure; plot(real(ev), imag(ev), 'x', 'MarkerSize', 8); grid on
xline(0, 'r--'); xlabel('Re'); ylabel('Im'); title('Two-stage, \alpha = 1, T_m = 20 ms');
