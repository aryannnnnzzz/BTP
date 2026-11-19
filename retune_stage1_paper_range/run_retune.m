function run_retune()
%RUN_RETUNE  Stage-1 retune within the base paper's gain range (separate study).
%  Run from THIS folder:   matlab -batch "run_retune"
%  Self-contained: uses only code\ in this folder and writes only into this
%  folder (stage1_result.mat, stage2_result.mat, CSVs, results_<date>\).
%  Nothing in the parent project is read or modified, except that
%  stage2_result_run4.mat / stage1_result_run4.mat here are copies of run 4.
%
%  Differences from run 4 (parent project, results_20261119_1549):
%   - stage-1 search range: Kp, Ki <= 200 and Kvi, Dvi <= 5 (paper Eq. 32),
%     Kd <= 20 (was 1000 / 1000 / 20 / 20, Kd <= 50)
%   - stage 2 is NOT retuned: the run-4 stage-2 parameters are reused, and
%     their robustness gates are re-checked with the new stage 1
here = fileparts(mfilename('fullpath'));
cd(here);
addpath(fullfile(here, 'code'));
% safety: every function must come from code\ in this folder
fns = {'sim_mg','mg_linearize','tune_stage1','run_validation','compare_paper', ...
       'stability_frozen','cost_stage1','ctrl_presets','mg_scenarios','loop_margin'};
for k = 1:numel(fns)
    w = which(fns{k});
    assert(startsWith(w, fullfile(here, 'code')), 'run_retune: %s resolves outside this folder (%s)', fns{k}, w);
end

outDir = fullfile(here, ['results_' datestr(now, 'yyyymmdd_HHMM')]);
mkdir(outDir);
diary(fullfile(outDir, 'log.txt'));
cleanup = onCleanup(@() diary('off'));
set(0, 'DefaultFigureVisible', 'off');

copyfile('stage2_result_run4.mat', 'stage2_result.mat');   % run-4 stage 2, unchanged

steps = {'tune_stage1', 'check_stage2_gates', 'run_validation', 'compare_paper', 'stability_frozen'};
for i = 1:numel(steps)
    fprintf('\n\n######## %s ########\n', steps{i});
    t0 = tic;
    if strcmp(steps{i}, 'check_stage2_gates')
        check_stage2_gates();
    else
        run_step(steps{i});
    end
    fprintf('######## %s done in %.1f min ########\n', steps{i}, toc(t0)/60);
    figs = findall(0, 'Type', 'figure');
    for k = 1:numel(figs)
        nm = get(figs(k), 'Name');
        if isempty(nm), nm = sprintf('fig%d', figs(k).Number); end
        nm = regexprep(nm, '[^\w]', '_');
        saveas(figs(k), fullfile(outDir, sprintf('%s_%s.png', steps{i}, nm)));
    end
    close all
end

copyfile('stage1_result.mat', outDir);
copyfile('stage2_result.mat', outDir);
copyfile('validation_results.csv', outDir);
copyfile('paper_comparison.csv', outDir);
set(0, 'DefaultFigureVisible', 'on');
fprintf('\nAll done. Results in %s\n', outDir);
end

function run_step(name)
run(name);
end

function check_stage2_gates()
% the run-4 stage 2 must still pass its own robustness gates with the new
% stage 1: zeta >= 0.10 at Tm = 20 ms and zeta >= 0.05 at Tm = 50 ms over
% alpha in {a_min, 1, a_hi} and kh in {kh_min, kh_max}
P = mg_params();  opt = sim_options();
load('stage1_result.mat', 'C1');  load('stage2_result.mat', 'S2');
ok = true;
fprintf('Run-4 stage 2 with the new stage 1 (alpha, kh, Tm -> max Re, zeta):\n');
for gate = [0.02 0.10; 0.05 0.05]'
    for a = unique([S2.a_min, 1, S2.a_hi])
        for kh = [S2.kh_min, S2.kh_max]
            [z, mr] = loop_margin(P, C1, a, kh, gate(1), opt);
            pass = mr < -1e-4 && z >= gate(2);
            ok = ok && pass;
            fprintf('  alpha %.2f  kh %.3f  Tm %2.0f ms -> Re %+.4g  zeta %.3f  (need >= %.2f)  %s\n', ...
                    a, kh, 1e3*gate(1), mr, z, gate(2), ternary(pass, 'ok', 'FAIL'));
        end
    end
end
if ok
    fprintf('All stage-2 gates pass: reusing the run-4 stage 2 is valid.\n');
else
    fprintf('WARNING: some stage-2 gates fail; stage 2 should be retuned for this stage 1.\n');
end
end

function s = ternary(c, a, b)
if c, s = a; else, s = b; end
end
