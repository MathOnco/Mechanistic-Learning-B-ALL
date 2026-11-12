function files = score_winner_subset_outer()
%SCORE_WINNER_SUBSET_OUTER Outer-test scores of the Delta-BIC winner alone.
%
% The winner is the most common exact subset across outer folds. Ties go to
% the higher mean inner-CV score, then smaller k, then smaller training BIC,
% then smaller mask. That one subset is refit on every outer-training fold.
% Balanced accuracy and AUC use only those held-out scores.

thisDir = fileparts(mfilename('fullpath'));
addpath(thisDir);
cfg = exhaustive_config();
cfg.outputDir = fullfile(thisDir, 'outputs', 'delta_bic_selection');
if ~isfolder(cfg.outputDir), mkdir(cfg.outputDir); end

selectionRoot = ['/Users/4477547/Moffitt Cancer Center Dropbox/' ...
    'Sadegh Marzban/Markov-Bcell-ALL/MechLearn Paper/Paper/' ...
    'Reviewer_Exhaustive_Nested_5Fold'];
sources = {
    "markov", fullfile(selectionRoot, 'outputs', 'delta_bic_selection')
    "flow", fullfile(selectionRoot, 'outputs_10fold', 'delta_bic_selection')
    };

files = strings(size(sources, 1), 1);
for v = 1:size(sources, 1)
    view = sources{v, 1};
    files(v) = scoreView(cfg, view, sources{v, 2});
    fprintf('Wrote %s\n', files(v));
end
end

function outFile = scoreView(cfg, view, selectionDir)
models = ["ridge","lasso","elasticnet"];
outcomes = ["BCR","MRDflow"];
summary = table();
for o = 1:numel(outcomes)
    C = exhaustive.loadCohort(cfg, outcomes(o), view);
    for m = 1:numel(models)
        selectionFile = fullfile(selectionDir, char(lower(string(outcomes(o))) + ...
            "__" + view + "__" + models(m) + "__deltaBIC2__fold_selections.csv"));
        selected = readtable(selectionFile, 'TextType', 'string', ...
            'VariableNamingRule', 'preserve');
        winner = modalSubset(selected);
        scores = refitWinner(cfg, C, outcomes(o), models(m), ...
            selected.fixedLambda(1), winner.mask, selected);
        point = exhaustive.metrics(C.y, scores);
        [ciLower, ciUpper] = bootstrapMetrics(C.y, scores, 2000, ...
            1729 + 100 * o + m);
        summary = [summary; table( ...
            view, outcomes(o), models(m), height(selected), ...
            string(winner.features), winner.count, ...
            point.balancedAccuracy, ciLower(1), ciUpper(1), ...
            point.auc, ciLower(2), ciUpper(2), ...
            'VariableNames', {'view','outcome','model','numberOuterSelections', ...
            'recommendedFinalSubset','outerSelectionsOfRecommendedSubset', ...
            'nestedBalancedAccuracy','balancedAccuracyCILower', ...
            'balancedAccuracyCIUpper','nestedAUC','aucCILower','aucCIUpper'})]; %#ok<AGROW>
        fprintf('%s %s %s: %s (%d/%d)  BA %.3f  AUC %.3f\n', ...
            view, outcomes(o), models(m), winner.features, winner.count, ...
            height(selected), point.balancedAccuracy, point.auc);
    end
end
outFile = string(fullfile(cfg.outputDir, ...
    char(view + "__deltaBIC2__final_subset_summary.csv")));
writetable(summary, outFile);
end

function winner = modalSubset(selected)
masks = unique(selected.selectedMask);
candidates = table();
for i = 1:numel(masks)
    hit = selected.selectedMask == masks(i);
    first = find(hit, 1);
    candidates(i, :) = table(masks(i), selected.k(first), ...
        selected.selectedFeatures(first), sum(hit), ...
        mean(selected.innerSelectionValue(hit)), ...
        mean(selected.trainingBIC(hit)), ...
        'VariableNames', {'mask','k','features','count', ...
        'meanInnerSelectionValue','meanTrainingBIC'});
end
candidates = sortrows(candidates, ...
    {'count','meanInnerSelectionValue','k','meanTrainingBIC','mask'}, ...
    {'descend','descend','ascend','ascend','ascend'});
winner = candidates(1, :);
end

function scores = refitWinner(cfg, C, outcome, model, lambda, mask, selected)
p = size(C.X, 2);
cols = exhaustive.indicesFromMask(mask, p);
assert(~isempty(cols), 'Winner mask %d has no features.', mask);
outcomeOffset = sum(double(char(outcome)));
outerSeed = cfg.randomSeed + outcomeOffset + 10000;
outer = exhaustive.stratifiedPartition(C.y, cfg.outerFolds, outerSeed);
assert(outer.NumTestSets == height(selected), ...
    'Outer-fold count does not match the saved Delta-BIC selections.');
scores = nan(numel(C.y), 1);
for f = 1:outer.NumTestSets
    tr = training(outer, f);
    te = test(outer, f);
    savedN = selected.nTest(selected.outerFold == f);
    assert(isscalar(savedN) && savedN == sum(te), ...
        'Rebuilt outer fold %d does not match the saved test count.', f);
    [Xtrain, Xtest] = exhaustive.prepareData(C.X(tr, :), C.X(te, :));
    scores(te) = exhaustive.fitPrepared(model, lambda, ...
        Xtrain(:, cols), C.y(tr), Xtest(:, cols), cfg);
end
assert(all(isfinite(scores)), 'Winner refit left a patient without a score.');
end

function [lower, upper] = bootstrapMetrics(y, score, iterations, seed)
y = double(y(:));
score = double(score(:));
positive = find(y == 1);
negative = find(y == 0);
rng(seed, 'twister');
boot = nan(iterations, 2);
for b = 1:iterations
    idx = [positive(randi(numel(positive), numel(positive), 1)); ...
           negative(randi(numel(negative), numel(negative), 1))];
    M = exhaustive.metrics(y(idx), score(idx));
    boot(b, :) = [M.balancedAccuracy, M.auc];
end
lower = prctile(boot, 2.5, 1);
upper = prctile(boot, 97.5, 1);
end
