%% COMPARE_PAPER  Like-for-like comparison with the base paper, Scenarios I-IV.
%  Sajadinia, J. Energy Storage 84 (2024) 110791, Tables 6, 8, 10, 12.
%  Same metrics as the paper (max overshoot, max undershoot, IAE) plus ITAE
%  (the paper's cost function, Eq. 21). All controllers run in the same model
%  with the exact VI block (Kvi*s + Dvi)/(1 + s*Tinv), dt = 0.2 ms.
%
%  Columns: paper-reported EO-FOPI and OSTIT2F-FOPI, then this model's
%  EO-FOPI (paper gains), OSTIT2F-FOPI (paper Table 5) and the two-stage
%  controller. "% vs EO" compares each controller with the EO-FOPI of the
%  SAME source, which is the fair comparison when the two models differ.
clear; clc; close all;
P = mg_params();
load('stage1_result.mat', 'C1');  load('stage2_result.mat', 'S2');
o = sim_options();  o.dt = 2e-4;  o.decim = 50;

Cp = ctrl_presets('eo_fopi');  Cp.Tf = 0;
Co = ctrl_presets('ostit2f');  Co.Tf = 0;

% paper values [OS US IAE]: rows EO-FOPI, OSTIT2F-FOPI
% (Tables 6, 10, 12 are shifted by one row in the PDF text; aligned with Table 8)
paper.I   = [0.00100 -0.00185 0.00722; 0.00077 -0.00133 0.00574];
paper.II  = [0.00069 -0.00154 0.00270; 0.00050 -0.00111 0.00213];
paper.III = [0.00043 -0.00239 0.00285; 0.00031 -0.00173 0.00225];
paper.IV  = [0.00076 -0.00162 0.00352; 0.00055 -0.00116 0.00280];

scen = {'I','II','III','IV'};
src  = {'Paper', 'Paper', 'This model', 'This model', 'This model'};
ctl  = {'EO-FOPI', 'OSTIT2F-FOPI', 'EO-FOPI', 'OSTIT2F-FOPI', 'Two-stage (proposed)'};
rows = {};
imp  = zeros(numel(scen), 3, 3);     % scenario x [OSTIT2F paper, OSTIT2F model, two-stage] x [OS US IAE]
for s = 1:numel(scen)
    sc = mg_scenarios(scen{s});
    rp = sim_mg(P, Cp, [], sc, o);
    ro = sim_mg(P, Co, [], sc, o);
    rd = sim_mg(P, C1, S2, sc, o);
    M  = [paper.(scen{s});
          rp.J.OS rp.J.US rp.J.IAE;
          ro.J.OS ro.J.US ro.J.IAE;
          rd.J.OS rd.J.US rd.J.IAE];
    itae = [NaN NaN rp.J.ITAE ro.J.ITAE rd.J.ITAE];
    ref  = [1 1 3 3 3];               % EO-FOPI row of the same source
    for k = 1:5
        pct = 100*(abs(M(k,:)) - abs(M(ref(k),:))) ./ abs(M(ref(k),:));
        rows(end+1,:) = {scen{s}, src{k}, ctl{k}, M(k,1), M(k,2), M(k,3), itae(k), ...
                         pct(1), pct(2), pct(3)}; %#ok<SAGROW>
    end
    imp(s,1,:) = 100*(abs(M(2,:)) - abs(M(1,:))) ./ abs(M(1,:));
    imp(s,2,:) = 100*(abs(M(4,:)) - abs(M(3,:))) ./ abs(M(3,:));
    imp(s,3,:) = 100*(abs(M(5,:)) - abs(M(3,:))) ./ abs(M(3,:));
end
T = cell2table(rows, 'VariableNames', {'Scenario','Source','Controller', ...
    'OS','US','IAE','ITAE','OS_pct_vs_EO','US_pct_vs_EO','IAE_pct_vs_EO'});
fprintf('Like-for-like comparison with the base paper (dt = 0.2 ms, exact VI block)\n');
disp(T);  writetable(T, 'paper_comparison.csv');

% this model's EO-FOPI against the paper's EO-FOPI (model agreement)
fprintf('\nModel agreement, this model''s EO-FOPI vs the paper''s EO-FOPI (%% difference):\n');
for s = 1:numel(scen)
    a = T{strcmp(T.Scenario, scen{s}) & strcmp(T.Source, 'This model') & strcmp(T.Controller, 'EO-FOPI'), {'OS','US','IAE'}};
    b = paper.(scen{s})(1,:);
    fprintf('  %-4s |OS| %+6.1f%%   |US| %+6.1f%%   IAE %+6.1f%%\n', scen{s}, 100*(abs(a) - abs(b))./abs(b));
end

% direct df -> ESS gain (high-frequency loop gain) of each controller
G = @(C) C.Kvi*(C.Kp + C.Kd*o.wh^C.mu)/P.Tinv;
fprintf('\nDirect df -> ESS power gain G (p.u./Hz): EO-FOPI %.1f, OSTIT2F-FOPI %.1f (at rest), two-stage %.1f\n', ...
        G(Cp), G(Co), G(C1));

% ---------- figure: improvement over each source's EO-FOPI ----------
figure('Name', 'paper_comparison');
lab = {'OS', 'US', 'IAE'};
for m = 1:3
    subplot(1,3,m);
    bar(categorical(scen), squeeze(imp(:,:,m)));  grid on
    ylabel('% change vs EO-FOPI (negative = better)');  title(lab{m});
    if m == 1
        legend({'OSTIT2F (paper-reported)', 'OSTIT2F (this model)', 'Two-stage (this model)'}, ...
               'Location', 'southoutside');
    end
end
