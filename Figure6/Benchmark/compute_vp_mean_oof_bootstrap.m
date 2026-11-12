function files = compute_vp_mean_oof_bootstrap(mode, models)
%COMPUTE_VP_MEAN_OOF_BOOTSTRAP Mean outer-OOF scores and bootstrap 95% CIs.
%
% The published VP selection checkpoints store only repeat-level metrics, not
% the held-out probabilities. This script re-runs the same nested CV design
% used for R3 VP feature stability, accumulates outer-OOF scores across
% repeats, then applies the same stratified patient/row bootstrap used for
% the original-data classifier benchmark.
%
%   files = compute_vp_mean_oof_bootstrap
%   files = compute_vp_mean_oof_bootstrap("smoke")
%   files = compute_vp_mean_oof_bootstrap("full")
%
% Progress is checkpointed after every repeat and can be resumed.

if nargin < 1, mode = "full"; end
if nargin < 2 || isempty(models)
    models = ["lasso","elasticnet"];
end
mode = lower(string(mode));
models = lower(string(models(:)'));
assert(ismember(mode, ["full","smoke"]), 'Mode must be "full" or "smoke".');
assert(all(ismember(models, ["ridge","lasso","elasticnet"])), ...
    'Models must be ridge, lasso, and/or elasticnet.');

thisDir = fileparts(mfilename('fullpath'));
figureRoot = fileparts(thisDir);
completeVPFile = fullfile(figureRoot, 'Data', 'PB_BM_merged_with_MRDflow_BCR.csv');
assert(isfile(completeVPFile), 'Missing VP table: %s', completeVPFile);
addpath(figureRoot);
addpath(fullfile(figureRoot, 'Functions'));

if mode == "smoke"
    outputDir = fullfile(figureRoot, 'outputs', 'smoke_vp_mean_oof_bootstrap');
    nRepeats = 1;
else
    outputDir = fullfile(figureRoot, 'outputs', 'specimen_folds');
    nRepeats = 100;
end
if ~isfolder(outputDir), mkdir(outputDir); end

cfg.nRepeats = nRepeats;
cfg.outerFolds = 5;
cfg.innerFolds = 5;
cfg.randomSeed = 1729;
cfg.lambdaGrid = logspace(-4, 2, 13);
cfg.elasticNetAlpha = 0.5;
cfg.primaryTuneMetric = "auc";
cfg.models = models;
cfg.vectorizedLambda = true;
cfg.bootstrapIterations = 2000;
warning('off', 'stats:lasso:MaxIterReached');
warning('off', 'stats:glmfit:IterationLimit');
if isequal(models, ["lasso","elasticnet"])
    fileTag = "";
else
    fileTag = "_" + strjoin(models, "_");
end

[Xall, featureNames, mrd, bcr, block] = loadCompleteVP(completeVPFile);
cohorts = struct( ...
    'name', {"BCR", "MRDflow"}, ...
    'eligible', {ismember(bcr, ["Negative","BCR::ABL1"]), ...
                 ismember(mrd, ["Negative","Positive (< 10^-3)", ...
                 "Positive (> 10^-3)"])}, ...
    'y', {double(bcr == "BCR::ABL1"), ...
          double(mrd == "Positive (> 10^-3)")});

files = struct();
for i = 1:numel(cohorts)
    keep = cohorts(i).eligible;
    X = Xall(keep, :);
    y = cohorts(i).y(keep);
    group = block(keep);
    label = "VP__" + cohorts(i).name + "__markov";
    stem = "vp_" + lower(cohorts(i).name) + "__markov" + fileTag;
    checkpointFile = fullfile(outputDir, char(stem + "__mean_oof_checkpoint.mat"));
    fprintf('\n%s: n=%d (+%d/-%d)\n', label, numel(y), sum(y==1), sum(y==0));

    signature = makeSignature(completeVPFile, label, X, y, featureNames, cfg);
    [meanOOF, modelNames] = accumulateMeanOOF( ...
        X, y, group, cfg, label, checkpointFile, signature);

    analysisLabel = label;
    summary = bootstrapSummary(y, meanOOF, modelNames, ...
        cfg.bootstrapIterations, cfg.randomSeed + 901, analysisLabel);
    bootFile = fullfile(outputDir, char(stem + "_bootstrap_95CI.csv"));
    oofFile = fullfile(outputDir, char(stem + "_mean_OOF_scores.csv"));
    writetable(summary, bootFile);
    oof = table((1:numel(y))', y, 'VariableNames', {'row','outcome'});
    for m = 1:numel(modelNames)
        oof.(matlab.lang.makeValidName(modelNames(m) + "_meanOOF")) = meanOOF(:, m);
    end
    writetable(oof, oofFile);
    files.(char(cohorts(i).name + "_bootstrap")) = string(bootFile);
    files.(char(cohorts(i).name + "_meanOOF")) = string(oofFile);
    fprintf('Wrote %s\n', bootFile);
end
end

function [meanOOF, modelNames] = accumulateMeanOOF( ...
        X, y, group, cfg, label, checkpointFile, signature)
n = size(X, 1);
nModels = numel(cfg.models);
modelNames = cfg.models(:);

state.completedRepeats = 0;
state.oofSum = zeros(n, nModels);
state.oofCount = zeros(n, nModels);
state.signature = signature;

if isfile(checkpointFile)
    saved = load(checkpointFile, 'state');
    assert(isfield(saved, 'state') && isequaln(saved.state.signature, signature), ...
        'Checkpoint signature mismatch: %s', checkpointFile);
    state = saved.state;
    fprintf('Resuming mean-OOF accumulation from repeat %d/%d.\n', ...
        state.completedRepeats, cfg.nRepeats);
end

y = double(y(:));
pending = (state.completedRepeats + 1):cfg.nRepeats;
if ~isempty(pending)
    functionsDir = fileparts(which('nestedCV'));
    if license('test', 'Distrib_Computing_Toolbox')
        pool = gcp('nocreate');
        if isempty(pool)
            pool = parpool;
        end
        pctRunOnAll(['addpath(''' strrep(functionsDir, '''', '''''') ''')']);
    end
    partial = cell(numel(pending), 1);
    parfor ii = 1:numel(pending)
        partial{ii} = oneSpecimenRepeat(X, y, group, cfg, pending(ii));
    end
    for ii = 1:numel(partial)
        state.oofSum = state.oofSum + partial{ii};
        state.oofCount = state.oofCount + 1;
    end
    state.completedRepeats = cfg.nRepeats;
    save(checkpointFile, 'state', '-v7.3');
    fprintf('%s mean-OOF repeats %d-%d complete.\n', ...
        label, pending(1), pending(end));
end

assert(all(state.oofCount(:) == cfg.nRepeats), ...
    'Incomplete OOF accumulation for %s.', label);
meanOOF = state.oofSum ./ state.oofCount;
end

function oofScores = oneSpecimenRepeat(X, y, group, cfg, repeatIndex)
warning('off', 'stats:lasso:MaxIterReached');
warning('off', 'stats:glmfit:IterationLimit');
n = numel(y);
nModels = numel(cfg.models);
outer = groupedStratifiedPartition(y, group, cfg.outerFolds, ...
    cfg.randomSeed + 10000 * repeatIndex);
oofScores = nan(n, nModels);
for f = 1:outer.NumTestSets
    tr = training(outer, f);
    te = test(outer, f);
    assertDisjointGroups(group(tr), group(te));
    for m = 1:nModels
        innerSeed = cfg.randomSeed + 100000 * repeatIndex + 100 * f + m;
        hp = tuneHyperparameter(cfg.models(m), X(tr, :), ...
            y(tr), cfg, innerSeed, group(tr));
        score = fitPredict(cfg.models(m), hp, ...
            X(tr, :), y(tr), X(te, :), cfg);
        oofScores(te, m) = score(:);
    end
end
assert(all(isfinite(oofScores), 'all'), ...
    'Nonfinite OOF scores on repeat %d.', repeatIndex);
end

function [X, featureNames, mrd, bcr, block] = loadCompleteVP(filename)
T = readtable(filename, 'VariableNamingRule', 'preserve', 'TextType', 'string');
featureNames = ["M11","M12","M13","M14","M21","M22","M23","M24", ...
    "M31","M32","M33","M34","M41","M42","M43","M44"];
required = ["ID","VirtualPatient",featureNames,"MRD(flow)","BCR::ABL1/like"];
missing = setdiff(required, string(T.Properties.VariableNames));
assert(isempty(missing), 'Complete VP table is missing: %s', strjoin(missing, ', '));
X = double(T{:, cellstr(featureNames)});
assert(all(isfinite(X), 'all'), 'VP Markov matrix has nonfinite values.');
mrd = strtrim(string(T{:, 'MRD(flow)'}));
bcr = strtrim(string(T{:, 'BCR::ABL1/like'}));
block = cumsum(double(T.VirtualPatient) == 1);
end

function signature = makeSignature(filename, label, X, y, featureNames, cfg)
signature.filename = string(filename);
signature.label = string(label);
signature.purpose = "mean_oof_bootstrap";
signature.n = size(X, 1);
signature.p = size(X, 2);
signature.nPositive = sum(y == 1);
signature.featureNames = featureNames;
signature.outerFolds = cfg.outerFolds;
signature.innerFolds = cfg.innerFolds;
signature.randomSeed = cfg.randomSeed;
signature.lambdaGrid = cfg.lambdaGrid;
signature.elasticNetAlpha = cfg.elasticNetAlpha;
signature.models = cfg.models;
signature.nRepeats = cfg.nRepeats;
signature.splitUnit = "source_specimen";
end
