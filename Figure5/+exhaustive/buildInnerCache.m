function cache = buildInnerCache(X, y, nFolds, seed, group)
%BUILDINNERCACHE Precompute training-only transforms for each inner fold.
%
% An optional group vector keeps every row with the same specimen ID in one
% inner fold. Omit it for the original row-level partition.

if nargin < 5 || isempty(group)
    cvp = exhaustive.stratifiedPartition(y, nFolds, seed);
else
    cvp = exhaustive.groupedStratifiedPartition(y, group, nFolds, seed);
end
cache.n = numel(y);
cache.y = double(y(:));
cache.folds = cell(cvp.NumTestSets, 1);
for f = 1:cvp.NumTestSets
    tr = training(cvp, f);
    te = test(cvp, f);
    if nargin >= 5 && ~isempty(group)
        exhaustive.assertDisjointGroups(group(tr), group(te));
    end
    [Xtr, Xte] = exhaustive.prepareData(X(tr, :), X(te, :));
    cache.folds{f} = struct( ...
        'trainIndex', find(tr), ...
        'testIndex', find(te), ...
        'Xtrain', Xtr, ...
        'Xtest', Xte, ...
        'ytrain', double(y(tr)), ...
        'ytest', double(y(te)));
end
end
