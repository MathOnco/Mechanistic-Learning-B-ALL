function auc = fastAUC(y, score)
%FASTAUC Exact rank-sum AUC with half credit for tied scores.
%
% This avoids the curve-construction overhead of perfcurve inside millions
% of subset evaluations while returning the same empirical ROC AUC.

y = double(y(:));
score = double(score(:));
nPos = sum(y == 1);
nNeg = sum(y == 0);
if nPos == 0 || nNeg == 0
    auc = NaN;
    return;
end

[sortedScore, order] = sort(score, 'ascend');
n = numel(score);
sortedRanks = (1:n)';
i = 1;
while i <= n
    j = i;
    while j < n && sortedScore(j + 1) == sortedScore(i)
        j = j + 1;
    end
    if j > i
        sortedRanks(i:j) = 0.5 * (i + j);
    end
    i = j + 1;
end
ranks = zeros(n, 1);
ranks(order) = sortedRanks;
auc = (sum(ranks(y == 1)) - nPos * (nPos + 1) / 2) / (nPos * nNeg);
auc = min(max(auc, 0), 1);
end
