function files = collect_vp_mrd_mean_lambdas()
%COLLECT_VP_MRD_MEAN_LAMBDAS Mean penalties from the specimen-level VP nested CV.
%
% The earlier out-of-fold run tuned lambda inside every outer training fold
% but did not save the chosen values. This repeats that tuning for flow MRD
% on all 100 virtual-patient rows per specimen. Lasso uses model index 1 and
% elastic net index 2, matching their joint run. Ridge uses model index 1,
% matching its separate run. The reported penalty is the mean of 500 choices.

thisDir = fileparts(mfilename('fullpath'));
outputDir = fullfile(thisDir, 'outputs');
if ~isfolder(outputDir), mkdir(outputDir); end
completeVPFile = fullfile(thisDir, 'Data', 'PB_BM_merged_with_MRDflow_BCR.csv');
assert(isfile(completeVPFile), 'Missing VP table: %s', completeVPFile);
addpath(thisDir);
addpath(fullfile(thisDir, 'Functions'));

cfg.nRepeats = 100;
cfg.outerFolds = 5;
cfg.innerFolds = 5;
cfg.randomSeed = 1729;
cfg.lambdaGrid = logspace(-4, 2, 13);
cfg.elasticNetAlpha = 0.5;
cfg.primaryTuneMetric = "auc";
cfg.vectorizedLambda = true;
cfg.models = ["lasso","elasticnet","ridge"];
% Inner-seed offsets used when these models were tuned for the OOF figures.
cfg.tuningModelIndex = [1, 2, 1];
warning('off', 'stats:lasso:MaxIterReached');
warning('off', 'stats:glmfit:IterationLimit');

[X, y, group] = loadMrdVP(completeVPFile);
fprintf('Flow MRD VP for lambda recovery: n=%d (+%d/-%d), %d specimens.\n', ...
    numel(y), sum(y == 1), sum(y == 0), numel(unique(group)));

checkpointFile = fullfile(outputDir, 'vp_mrd_selected_lambdas_checkpoint.mat');
signature = struct('purpose', "vp_mrd_mean_lambda", 'n', numel(y), ...
    'nRepeats', cfg.nRepeats, 'outerFolds', cfg.outerFolds, ...
    'innerFolds', cfg.innerFolds, 'randomSeed', cfg.randomSeed, ...
    'models', cfg.models, 'tuningModelIndex', cfg.tuningModelIndex, ...
    'splitUnit', "source_specimen");
if isfile(checkpointFile)
    saved = load(checkpointFile, 'lambdas', 'signature');
    assert(isequaln(saved.signature, signature), ...
        'Lambda checkpoint does not match the current VP nested design.');
    lambdas = saved.lambdas;
    fprintf('Loaded saved lambda choices from %s\n', checkpointFile);
else
    functionsDir = fileparts(which('tuneHyperparameter'));
    pool = gcp('nocreate');
    if isempty(pool), pool = parpool; end
    pctRunOnAll(['addpath(''' strrep(functionsDir, '''', '''''') ''')']);
    fprintf('Parallel pool: %d workers\n', pool.NumWorkers);
    partial = cell(cfg.nRepeats, 1);
    parfor r = 1:cfg.nRepeats
        partial{r} = oneRepeatLambdas(X, y, group, cfg, r);
    end
    lambdas = cat(1, partial{:});
    save(checkpointFile, 'lambdas', 'signature', '-v7.3');
end

choiceFile = fullfile(outputDir, 'vp_mrd_selected_lambdas.csv');
choices = table();
for m = 1:numel(cfg.models)
    values = lambdas(:, m);
    repeat = repelem((1:cfg.nRepeats)', cfg.outerFolds);
    fold = repmat((1:cfg.outerFolds)', cfg.nRepeats, 1);
    block = table(repmat(cfg.models(m), numel(values), 1), repeat, fold, values, ...
        'VariableNames', {'model','repeat','outerFold','lambda'});
    choices = [choices; block]; %#ok<AGROW>
end
writetable(choices, choiceFile);

summary = table();
for m = 1:numel(cfg.models)
    values = lambdas(:, m);
    summary = [summary; table(cfg.models(m), numel(values), mean(values), ...
        median(values), std(values), ...
        'VariableNames', {'model','nOuterFits','meanLambda','medianLambda','sdLambda'})]; %#ok<AGROW>
end
summaryFile = fullfile(outputDir, 'vp_mrd_mean_lambdas.csv');
writetable(summary, summaryFile);
disp(summary);
fprintf('Wrote %s\n', summaryFile);

files.choices = string(choiceFile);
files.summary = string(summaryFile);
end

function chosen = oneRepeatLambdas(X, y, group, cfg, repeatIndex)
warning('off', 'stats:lasso:MaxIterReached');
warning('off', 'stats:glmfit:IterationLimit');
outer = groupedStratifiedPartition(y, group, cfg.outerFolds, ...
    cfg.randomSeed + 10000 * repeatIndex);
chosen = nan(outer.NumTestSets, numel(cfg.models));
for f = 1:outer.NumTestSets
    tr = training(outer, f);
    te = test(outer, f);
    assertDisjointGroups(group(tr), group(te));
    for m = 1:numel(cfg.models)
        innerSeed = cfg.randomSeed + 100000 * repeatIndex + 100 * f + ...
            cfg.tuningModelIndex(m);
        chosen(f, m) = tuneHyperparameter(cfg.models(m), ...
            X(tr, :), y(tr), cfg, innerSeed, group(tr));
    end
end
end

function [X, y, block] = loadMrdVP(filename)
T = readtable(filename, 'VariableNamingRule', 'preserve', 'TextType', 'string');
featureNames = ["M11","M12","M13","M14","M21","M22","M23","M24", ...
    "M31","M32","M33","M34","M41","M42","M43","M44"];
X = double(T{:, cellstr(featureNames)});
mrd = strtrim(string(T{:, 'MRD(flow)'}));
keep = ismember(mrd, ["Negative","Positive (< 10^-3)","Positive (> 10^-3)"]);
block = cumsum(double(T.VirtualPatient) == 1);
X = X(keep, :);
block = block(keep);
y = double(mrd(keep) == "Positive (> 10^-3)");
assert(all(isfinite(X), 'all'), 'VP Markov matrix has nonfinite values.');
end
