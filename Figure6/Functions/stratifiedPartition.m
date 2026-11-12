function cvp = stratifiedPartition(y, requestedFolds, seed)
%STRATIFIEDPARTITION Reproducible K-fold split with feasible class counts.

y = double(y(:));
nPos = sum(y == 1);
nNeg = sum(y == 0);
k = min([requestedFolds, nPos, nNeg]);
if k < 2
    error('mechlearn:InsufficientClassCount', ...
        'At least two observations from each class are required; found %d and %d.', nPos, nNeg);
end
rng(seed, 'twister');
Ycat = categorical(y, [0 1], {'Negative','Positive'});
cvp = cvpartition(Ycat, 'KFold', k);
end
