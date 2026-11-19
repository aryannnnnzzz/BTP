function run_all()
%RUN_ALL  Run the whole re-tune + validation pipeline without supervision.
%  Works from the MATLAB desktop or from a terminal (e.g. Claude Code):
%      matlab -batch "run_all"
%  Saves everything into a time-stamped folder results_<date>:
%      log.txt (all printed output), all figures as .png, the .mat results
%      and validation_results.csv.
%  Old stage1/stage2 results are backed up as *_v1.mat before re-tuning.

outDir = ['results_' datestr(now, 'yyyymmdd_HHMM')];
mkdir(outDir);
diary(fullfile(outDir, 'log.txt'));
cleanup = onCleanup(@() diary('off'));
set(0, 'DefaultFigureVisible', 'off');     % no windows needed in batch mode

files = {'stage1_result.mat', 'stage2_result.mat'};
for k = 1:numel(files)
    bak = strrep(files{k}, '.mat', '_v1.mat');
    if exist(files{k}, 'file') && ~exist(bak, 'file')
        copyfile(files{k}, bak);
    end
end

steps = {'tune_stage1', 'tune_stage2', 'run_validation', 'compare_paper', 'stability_frozen'};
for i = 1:numel(steps)
    fprintf('\n\n######## %s ########\n', steps{i});
    t0 = tic;
    run_step(steps{i});                    % own workspace: the scripts' "clear" is harmless
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
if exist('paper_comparison.csv', 'file')
    copyfile('paper_comparison.csv', outDir);
end
if exist('validation_results.csv', 'file')
    copyfile('validation_results.csv', outDir);
end
set(0, 'DefaultFigureVisible', 'on');
fprintf('\nAll done. Results in %s\n', outDir);
end

function run_step(name)
run(name);
end
