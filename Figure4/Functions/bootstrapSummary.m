function T = bootstrapSummary(y, scores, modelNames, nBoot, seed, analysisLabel)
%BOOTSTRAPSUMMARY Stratified-bootstrap 95% CIs for OOF performance.

y = double(y(:));
scores = double(scores);
if isvector(scores)
    scores = scores(:);
end
modelNames = string(modelNames(:));
pos = find(y == 1);
neg = find(y == 0);
nModels = size(scores, 2);
rng(seed, 'twister');

metricNames = ["accuracy","balancedAccuracy","sensitivity","specificity","auc"];
boot = nan(nBoot, numel(metricNames), nModels);
for b = 1:nBoot
    idx = [pos(randi(numel(pos), numel(pos), 1)); ...
           neg(randi(numel(neg), numel(neg), 1))];
    for m = 1:nModels
        M = metrics(y(idx), scores(idx, m));
        boot(b, :, m) = [M.accuracy, M.balancedAccuracy, ...
            M.sensitivity, M.specificity, M.auc];
    end
end

T = table();
row = 0;
for m = 1:nModels
    point = metrics(y, scores(:, m));
    pointValues = [point.accuracy, point.balancedAccuracy, ...
        point.sensitivity, point.specificity, point.auc];
    for j = 1:numel(metricNames)
        row = row + 1;
        ci = prctile(boot(:, j, m), [2.5 97.5]);
        T(row, :) = table(string(analysisLabel), modelNames(m), metricNames(j), ...
            pointValues(j), ci(1), ci(2), point.n, point.nPositive, point.nNegative, ...
            'VariableNames', {'analysis','model','metric','estimate','ciLower','ciUpper', ...
            'n','nPositive','nNegative'});
    end
end
end
