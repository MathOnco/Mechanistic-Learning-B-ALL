function run_exhaustive_nested_10fold(outcomes, models, views)
%RUN_EXHAUSTIVE_NESTED_10FOLD Run best-k exhaustive searches under nested CV.

clc;
thisDir = fileparts(mfilename('fullpath'));
addpath(thisDir);
cfg = exhaustive_config();
if nargin >= 1 && ~isempty(outcomes), cfg.outcomes = string(outcomes); end
if nargin >= 2 && ~isempty(models), cfg.models = lower(string(models)); end
if nargin >= 3 && ~isempty(views), cfg.views = lower(string(views)); end

fullSubsetCount = sum(arrayfun(@(k)nchoosek(16, k), 1:16));
assert(fullSubsetCount == 2^16 - 1 && fullSubsetCount == 65535, ...
    'The 16-transition subset count must equal 65,535.');
requestedSubsetCount = 0;
for v = 1:numel(cfg.views)
    pView = numberOfFeatures(cfg, cfg.views(v));
    kView = min(cfg.maxK, pView);
    requestedSubsetCount = requestedSubsetCount + ...
        sum(arrayfun(@(k)nchoosek(pView, k), 1:kView));
end
estimatedFits = numel(cfg.outcomes) * numel(cfg.models) * cfg.nRepeats * ...
    cfg.outerFolds * cfg.innerFolds * requestedSubsetCount;
% Parallel-pool startup costs more than the entire four-feature Flow search.
if requestedSubsetCount <= 100
    cfg.useParallel = false;
end

fprintf('Ten-fold nested exhaustive transition search\n');
fprintf('Output: %s\n', cfg.outputDir);
fprintf('Selection metric: %s\n', cfg.selectionMetric);
if cfg.selectionMetric == "combo"
    fprintf('Score weights (accuracy, balanced accuracy, AUC): %.2f, %.2f, %.2f\n', ...
        cfg.comboWeights(1), cfg.comboWeights(2), cfg.comboWeights(3));
end
if isfinite(cfg.maxK), fprintf('Maximum requested k: %d\n', cfg.maxK); end
fprintf('All non-empty subsets for 16 transitions: %s\n', ...
    char(string(fullSubsetCount)));
for v = 1:numel(cfg.views)
    pView = numberOfFeatures(cfg, cfg.views(v));
    fprintf('%s predictors: %d features, %d non-empty subsets\n', ...
        cfg.views(v), pView, 2^pView - 1);
end
fprintf('Approximate inner model fits for this invocation: %s\n', ...
    char(string(estimatedFits)));

if cfg.useParallel && license('test', 'Distrib_Computing_Toolbox')
    pool = gcp('nocreate');
    if isempty(pool)
        try
            pool = parpool;
        catch ME
            warning('exhaustive:ParallelUnavailable', ...
                'Could not start a parallel pool; continuing serially: %s', ME.message);
            cfg.useParallel = false;
        end
    end
    if ~isempty(pool)
        fprintf('Parallel pool: %d workers\n', pool.NumWorkers);
    end
else
    cfg.useParallel = false;
end

allBestFeatures = table();
for v = 1:numel(cfg.views)
    view = cfg.views(v);
    for o = 1:numel(cfg.outcomes)
        outcome = cfg.outcomes(o);
        C = exhaustive.loadCohort(cfg, outcome, view);
        fprintf('\n%s | %s: N=%d, positive=%d, negative=%d, predictors=%d\n', ...
            view, outcome, numel(C.y), sum(C.y == 1), sum(C.y == 0), ...
            size(C.X, 2));

        for m = 1:numel(cfg.models)
            model = cfg.models(m);
            lambda = exhaustive.lambdaFor(cfg, view, outcome, model);
            fprintf('\n  %s, fixed lambda %.15g\n', model, lambda);
            R = exhaustive.runEndpointModel(C, outcome, model, lambda, cfg);
            bestFeatures = exhaustive.saveResults(R, C, cfg);
            allBestFeatures = [allBestFeatures; bestFeatures]; %#ok<AGROW>
        end
    end
end

if isscalar(cfg.views)
    combinedName = char(cfg.views(1) + "__all_models_best_features_and_accuracies.csv");
else
    combinedName = 'all_views_models_best_features_and_accuracies.csv';
end
writetable(allBestFeatures, fullfile(cfg.outputDir, combinedName));

fprintf('\nCompleted. Results are in %s\n', cfg.outputDir);
end

function p = numberOfFeatures(cfg, view)
switch lower(string(view))
    case "markov"
        p = numel(cfg.markovFeatureNames);
    case "flow"
        p = numel(cfg.flowFeatureNames);
    otherwise
        error('exhaustive:UnknownView', 'Unknown predictor view: %s', view);
end
end
