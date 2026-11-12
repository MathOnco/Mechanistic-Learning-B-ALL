function R = nestedCV(X, y, featureNames, modelNames, cfg, analysisLabel, group)
%NESTEDCV Repeated nested stratified CV with stored OOF predictions.
%
% Pass group to keep every row with the same ID entirely in training or
% entirely in test, in both the outer and inner partitions. Omit group for
% the original row-level partition.

y = double(y(:));
featureNames = string(featureNames(:)');
modelNames = string(modelNames(:)');
if nargin < 7
    group = [];
end
group = double(group(:));
n = numel(y);
if ~isempty(group) && numel(group) ~= n
    error('mechlearn:GroupLengthMismatch', ...
        'Expected one group ID per observation; found %d IDs and %d rows.', ...
        numel(group), n);
end
nModels = numel(modelNames);
nFeatures = size(X, 2);

oofScores = nan(n, cfg.nRepeats, nModels);
repeatMetrics = table();
selectionCounts = zeros(nFeatures, nModels);
fitCounts = zeros(1, nModels);
hyperparameters = cell(cfg.nRepeats, nModels);

row = 0;
for r = 1:cfg.nRepeats
    outerSeed = cfg.randomSeed + 10000*r;
    if isempty(group)
        outer = stratifiedPartition(y, cfg.outerFolds, outerSeed);
    else
        outer = groupedStratifiedPartition( ...
            y, group, cfg.outerFolds, outerSeed);
    end
    hpRepeat = cell(outer.NumTestSets, nModels);
    for f = 1:outer.NumTestSets
        tr = training(outer, f);
        te = test(outer, f);
        if ~isempty(group)
            assertDisjointGroups(group(tr), group(te));
        end
        for m = 1:nModels
            modelName = modelNames(m);
            seed = cfg.randomSeed + 100000*r + 100*f + m;
            if isempty(group)
                hp = tuneHyperparameter(modelName, X(tr, :), y(tr), cfg, seed);
            else
                hp = tuneHyperparameter(modelName, X(tr, :), y(tr), ...
                    cfg, seed, group(tr));
            end
            [score, selected] = fitPredict(modelName, hp, ...
                X(tr, :), y(tr), X(te, :), cfg);
            oofScores(te, r, m) = score;
            selectionCounts(:, m) = selectionCounts(:, m) + selected(:);
            fitCounts(m) = fitCounts(m) + 1;
            hpRepeat{f, m} = hp;
        end
    end

    for m = 1:nModels
        M = metrics(y, oofScores(:, r, m));
        row = row + 1;
        repeatMetrics(row, :) = table(string(analysisLabel), r, modelNames(m), ...
            M.n, M.nPositive, M.nNegative, M.accuracy, M.balancedAccuracy, ...
            M.sensitivity, M.specificity, M.auc, ...
            'VariableNames', {'analysis','repeat','model','n','nPositive','nNegative', ...
            'accuracy','balancedAccuracy','sensitivity','specificity','auc'});
        hyperparameters{r, m} = cell2mat(hpRepeat(:, m));
    end
end

if any(isnan(oofScores), 'all')
    error('mechlearn:IncompleteOOF', 'At least one observation lacks an outer-fold prediction.');
end

selectionFrequency = selectionCounts ./ max(1, fitCounts);
selectionTable = table(featureNames', 'VariableNames', {'feature'});
for m = 1:nModels
    selectionTable.(matlab.lang.makeValidName(modelNames(m) + "_frequency")) = ...
        selectionFrequency(:, m);
end

R.analysisLabel = string(analysisLabel);
R.y = y;
R.featureNames = featureNames;
R.modelNames = modelNames;
R.oofScores = oofScores;
R.meanOOFScores = squeeze(mean(oofScores, 2));
if nModels == 1
    R.meanOOFScores = R.meanOOFScores(:);
end
R.repeatMetrics = repeatMetrics;
R.selectionTable = selectionTable;
R.hyperparameters = hyperparameters;
end
