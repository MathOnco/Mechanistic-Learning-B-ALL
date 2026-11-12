function M = metrics(y, score, threshold)
%METRICS Binary classification metrics from positive-class scores.

if nargin < 3
    threshold = 0.5;
end
y = double(y(:));
score = double(score(:));
valid = isfinite(y) & isfinite(score);
y = y(valid);
score = score(valid);
pred = score >= threshold;

M.n = numel(y);
M.nPositive = sum(y == 1);
M.nNegative = sum(y == 0);
M.accuracy = mean(pred == y);
M.sensitivity = sum(pred == 1 & y == 1) / max(1, sum(y == 1));
M.specificity = sum(pred == 0 & y == 0) / max(1, sum(y == 0));
M.balancedAccuracy = 0.5 * (M.sensitivity + M.specificity);
if numel(unique(y)) == 2
    [~, ~, ~, M.auc] = perfcurve(y, score, 1);
else
    M.auc = NaN;
end
end
