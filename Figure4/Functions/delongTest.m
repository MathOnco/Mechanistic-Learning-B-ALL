function [difference, pValue, covariance] = delongTest(y, scoreA, scoreB)
%DELONGTEST Paired nonparametric DeLong test for two ROC AUCs.
% Applied to each sample's mean repeated-OOF score as a secondary analysis.

y = double(y(:));
scores = [double(scoreA(:)), double(scoreB(:))];
posScores = scores(y == 1, :);
negScores = scores(y == 0, :);
m = size(posScores, 1);
n = size(negScores, 1);
if m < 2 || n < 2
    difference = NaN; pValue = NaN; covariance = nan(2); return;
end

v10 = zeros(m, 2);
v01 = zeros(n, 2);
for model = 1:2
    comparisons = double(posScores(:, model) > negScores(:, model)') + ...
        0.5 * double(posScores(:, model) == negScores(:, model)');
    v10(:, model) = mean(comparisons, 2);
    v01(:, model) = mean(comparisons, 1)';
end
aucs = mean(v10, 1);
covariance = cov(v10, 1) / m + cov(v01, 1) / n;
contrast = [1 -1];
variance = contrast * covariance * contrast';
difference = aucs(1) - aucs(2);
if variance <= 0 || ~isfinite(variance)
    pValue = NaN;
else
    z = difference / sqrt(variance);
    pValue = 2 * normcdf(-abs(z));
end
end
