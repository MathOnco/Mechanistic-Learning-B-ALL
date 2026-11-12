function [selectionValue, metricVector, failed] = evaluateSubsetInnerCV(cache, cols, model, lambda, cfg)
%EVALUATESUBSETINNERCV Return inner out-of-fold performance for one subset.

scores = nan(cache.n, 1);
failed = false;

try
    for f = 1:numel(cache.folds)
        F = cache.folds{f};
        foldScore = fitPrepared(model, lambda, ...
            F.Xtrain(:, cols), F.ytrain, F.Xtest(:, cols), cfg);
        scores(F.testIndex) = foldScore;
    end
    if any(~isfinite(scores))
        error('exhaustive:IncompleteInnerOOF', 'Incomplete inner-fold predictions.');
    end
    M = metrics(cache.y, scores);
    combo = cfg.comboWeights(1) * M.accuracy + ...
        cfg.comboWeights(2) * M.balancedAccuracy + ...
        cfg.comboWeights(3) * M.auc;
    metricVector = [M.accuracy, M.balancedAccuracy, M.auc, combo];
    switch lower(cfg.selectionMetric)
        case "accuracy"
            selectionValue = M.accuracy;
        case "balancedaccuracy"
            selectionValue = M.balancedAccuracy;
        case "auc"
            selectionValue = M.auc;
        case "combo"
            selectionValue = combo;
        otherwise
            error('exhaustive:UnknownMetric', ...
                'Unknown selection metric: %s', cfg.selectionMetric);
    end
catch
    failed = true;
    selectionValue = -Inf;
    metricVector = [NaN, NaN, NaN, NaN];
end
end
