function R = runEndpointModel(C, outcome, model, lambda, cfg)
%RUNENDPOINTMODEL Nested exhaustive best-k path for one endpoint and model.

outcome = string(outcome);
model = lower(string(model));
X = double(C.X);
y = double(C.y(:));
[n, p] = size(X);
K = min(cfg.maxK, p);
if p ~= numel(C.featureNames) || p < 1 || p > 16
    error('exhaustive:InvalidFeatureCount', ...
        'Expected 1 to 16 named predictors; found %d predictors and %d names.', ...
        p, numel(C.featureNames));
end
if K < 1 || K > p
    error('exhaustive:InvalidMaxK', 'maxK must be between 1 and %d.', p);
end

view = lower(string(C.view));
stem = lower(char(outcome + "__" + view + "__" + model));
selectionTag = sprintf('%s_w%03d_%03d_%03d', char(cfg.selectionMetric), ...
    round(100 * cfg.comboWeights(1)), round(100 * cfg.comboWeights(2)), ...
    round(100 * cfg.comboWeights(3)));
checkpointFile = fullfile(cfg.outputDir, ...
    [stem '__' selectionTag '__checkpoint.mat']);
meta = struct('view', view, 'outcome', outcome, 'model', model, 'lambda', lambda, ...
    'n', n, 'p', p, 'K', K, 'nRepeats', cfg.nRepeats, ...
    'outerFolds', cfg.outerFolds, 'innerFolds', cfg.innerFolds, ...
    'selectionMetric', cfg.selectionMetric, ...
    'comboWeights', cfg.comboWeights, 'randomSeed', cfg.randomSeed);

state = initializeState(n, cfg.nRepeats, cfg.outerFolds, K);
if cfg.resumeFromCheckpoint && isfile(checkpointFile)
    saved = load(checkpointFile, 'state', 'meta');
    validateCheckpoint(saved.meta, meta, checkpointFile);
    state = saved.state;
    fprintf('  Resuming checkpoint: %s\n', checkpointFile);
end

hasGroup = isfield(C, 'group') && ~isempty(C.group);
if hasGroup && numel(C.group) ~= n
    error('exhaustive:GroupLengthMismatch', ...
        'Expected one specimen ID per row; found %d IDs and %d rows.', ...
        numel(C.group), n);
end
meta.splitUnit = "row";
if hasGroup, meta.splitUnit = "source_specimen"; end

outcomeOffset = sum(double(char(outcome)));
for r = 1:cfg.nRepeats
    outerSeed = cfg.randomSeed + outcomeOffset + 10000 * r;
    if hasGroup
        outer = groupedStratifiedPartition( ...
            y, C.group, cfg.outerFolds, outerSeed);
    else
        outer = stratifiedPartition(y, cfg.outerFolds, outerSeed);
    end
    if outer.NumTestSets ~= cfg.outerFolds
        error('exhaustive:OuterFoldMismatch', ...
            'Requested %d outer folds but only %d were feasible.', ...
            cfg.outerFolds, outer.NumTestSets);
    end

    for f = 1:outer.NumTestSets
        tr = training(outer, f);
        te = test(outer, f);
        if hasGroup
            assertDisjointGroups(C.group(tr), C.group(te));
        end
        state.outerFoldIndex(te, r) = f;
        if all(state.completed(r, f, :), 'all')
            fprintf('  repeat %d, outer fold %d already complete\n', r, f);
            continue;
        end

        innerSeed = cfg.randomSeed + outcomeOffset + 100000 * r + 100 * f;
        if hasGroup
            innerCache = buildInnerCache(X(tr, :), y(tr), ...
                cfg.innerFolds, innerSeed, C.group(tr));
        else
            innerCache = buildInnerCache(X(tr, :), y(tr), ...
                cfg.innerFolds, innerSeed);
        end
        [XouterTrain, XouterTest] = prepareData(X(tr, :), X(te, :));

        for k = 1:K
            if state.completed(r, f, k), continue; end
            fprintf('  repeat %d/%d, outer fold %d/%d, k=%d: %d subsets\n', ...
                r, cfg.nRepeats, f, outer.NumTestSets, k, nchoosek(p, k));

            best = searchBestByK(innerCache, k, p, model, lambda, cfg);
            cols = best.cols;
            [testScore, activeLocal, modelInfo] = fitPrepared( ...
                model, lambda, XouterTrain(:, cols), y(tr), ...
                XouterTest(:, cols), cfg);
            trainInfo = informationCriteria(modelInfo, ...
                XouterTrain(:, cols), y(tr), k, p, cfg);

            activeCols = cols(activeLocal);
            state.oofScores(te, r, k) = testScore;
            state.selectedMasks(r, f, k) = maskFromIndices(cols);
            state.activeMasks(r, f, k) = maskFromIndices(activeCols);
            state.innerMetrics(r, f, k, :) = best.metricVector;
            state.innerSelectionValue(r, f, k) = best.selectionValue;
            state.nCombinations(r, f, k) = best.nCombinations;
            state.nFailed(r, f, k) = best.nFailed;
            state.elapsedSeconds(r, f, k) = best.elapsedSeconds;
            state.logLikelihood(r, f, k) = trainInfo.logLikelihood;
            state.df(r, f, k) = trainInfo.df;
            state.bic(r, f, k) = trainInfo.bic;
            state.ebicGamma05(r, f, k) = trainInfo.ebicGamma05;
            state.ebicGamma1(r, f, k) = trainInfo.ebicGamma1;
            state.completed(r, f, k) = true;

            fprintf('    best: %s | inner %s=%.4f | %.1f s | failures=%d\n', ...
                strjoin(C.featureNames(cols), ','), cfg.selectionMetric, ...
                best.selectionValue, best.elapsedSeconds, best.nFailed);
            if cfg.saveCheckpointEveryK
                save(checkpointFile, 'state', 'meta', '-v7.3');
            end
        end
    end
end

if any(~state.completed, 'all') || any(~isfinite(state.oofScores), 'all')
    error('exhaustive:IncompleteResult', ...
        'The nested search ended with incomplete selections or predictions.');
end
save(checkpointFile, 'state', 'meta', '-v7.3');

R = state;
R.meta = meta;
R.outcome = outcome;
R.view = view;
R.model = model;
R.lambda = lambda;
R.featureNames = C.featureNames;
R.meanOOFScores = squeeze(mean(state.oofScores, 2));
if K == 1, R.meanOOFScores = R.meanOOFScores(:); end
R.checkpointFile = string(checkpointFile);
end

function state = initializeState(n, nRepeats, nFolds, K)
state.oofScores = nan(n, nRepeats, K);
state.outerFoldIndex = zeros(n, nRepeats);
state.selectedMasks = zeros(nRepeats, nFolds, K, 'uint16');
state.activeMasks = zeros(nRepeats, nFolds, K, 'uint16');
state.innerMetrics = nan(nRepeats, nFolds, K, 4);
state.innerSelectionValue = nan(nRepeats, nFolds, K);
state.nCombinations = zeros(nRepeats, nFolds, K);
state.nFailed = zeros(nRepeats, nFolds, K);
state.elapsedSeconds = nan(nRepeats, nFolds, K);
state.logLikelihood = nan(nRepeats, nFolds, K);
state.df = nan(nRepeats, nFolds, K);
state.bic = nan(nRepeats, nFolds, K);
state.ebicGamma05 = nan(nRepeats, nFolds, K);
state.ebicGamma1 = nan(nRepeats, nFolds, K);
state.completed = false(nRepeats, nFolds, K);
end

function validateCheckpoint(old, current, filename)
fields = {'view','outcome','model','lambda','n','p','K','nRepeats', ...
    'outerFolds','innerFolds','selectionMetric','comboWeights','randomSeed'};
for i = 1:numel(fields)
    name = fields{i};
    if ~isfield(old, name) || ~isequal(old.(name), current.(name))
        error('exhaustive:CheckpointMismatch', ...
            'Checkpoint %s does not match current setting %s.', filename, name);
    end
end
end
