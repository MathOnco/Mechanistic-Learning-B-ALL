function best = tuneHyperparameter(modelName, X, y, cfg, seed, group)
%TUNEHYPERPARAMETER Select model hyperparameters using inner CV only.
%
% An optional group vector keeps every row with the same ID entirely in the
% inner training fold or entirely in the inner test fold.

modelName = lower(string(modelName));
switch modelName
    case {"ridge", "lasso", "elasticnet"}
        candidates = cfg.lambdaGrid;
    case "randomforest"
        candidates = cfg.rfMinLeafGrid;
    otherwise
        error('mechlearn:UnknownModel', 'Unknown classifier: %s', modelName);
end

if nargin < 6 || isempty(group)
    cvp = stratifiedPartition(y, cfg.innerFolds, seed);
else
    cvp = groupedStratifiedPartition(y, group, cfg.innerFolds, seed);
end
useLambdaPath = isfield(cfg, 'vectorizedLambda') && cfg.vectorizedLambda && ...
    ismember(modelName, ["lasso","elasticnet"]) && numel(candidates) > 1;
candidateScores = -inf(size(candidates));
if useLambdaPath
    oof = nan(numel(y), numel(candidates));
    for f = 1:cvp.NumTestSets
        tr = training(cvp, f);
        te = test(cvp, f);
        oof(te, :) = fitPredict(modelName, candidates(:)', ...
            X(tr, :), y(tr), X(te, :), cfg);
    end
    for c = 1:numel(candidates)
        candidateScores(c) = selectionScore(y, oof(:, c), cfg);
    end
else
    for c = 1:numel(candidates)
        oof = nan(numel(y), 1);
        for f = 1:cvp.NumTestSets
            tr = training(cvp, f);
            te = test(cvp, f);
            oof(te) = fitPredict(modelName, candidates(c), ...
                X(tr, :), y(tr), X(te, :), cfg);
        end
        candidateScores(c) = selectionScore(y, oof, cfg);
    end
end

maxScore = max(candidateScores);
ties = find(abs(candidateScores - maxScore) < 1e-12);
% Prefer the simpler candidate on an exact tie: a larger leaf for random
% forests or stronger regularization for a penalized linear model.
[~, ii] = max(candidates(ties));
best = candidates(ties(ii));
end

function value = selectionScore(y, score, cfg)
M = metrics(y, score);
if cfg.primaryTuneMetric == "balancedaccuracy"
    value = M.balancedAccuracy;
else
    value = M.auc;
end
end
