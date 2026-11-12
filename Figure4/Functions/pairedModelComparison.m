function T = pairedModelComparison(y, scoreA, scoreB, labelA, labelB, nBoot, seed, analysisLabel)
%PAIREDMODELCOMPARISON Paired bootstrap differences and secondary DeLong test.

y = double(y(:));
scoreA = double(scoreA(:));
scoreB = double(scoreB(:));
pos = find(y == 1);
neg = find(y == 0);
rng(seed, 'twister');

metricNames = ["accuracy","balancedAccuracy","auc"];
pointA = metrics(y, scoreA);
pointB = metrics(y, scoreB);
valuesA = [pointA.accuracy, pointA.balancedAccuracy, pointA.auc];
valuesB = [pointB.accuracy, pointB.balancedAccuracy, pointB.auc];
diffBoot = nan(nBoot, numel(metricNames));

for b = 1:nBoot
    idx = [pos(randi(numel(pos), numel(pos), 1)); ...
           neg(randi(numel(neg), numel(neg), 1))];
    A = metrics(y(idx), scoreA(idx));
    B = metrics(y(idx), scoreB(idx));
    diffBoot(b, :) = [A.accuracy-B.accuracy, ...
        A.balancedAccuracy-B.balancedAccuracy, A.auc-B.auc];
end

T = table();
for j = 1:numel(metricNames)
    ci = prctile(diffBoot(:, j), [2.5 97.5]);
    pBoot = 2 * min(mean(diffBoot(:, j) <= 0), mean(diffBoot(:, j) >= 0));
    T(j, :) = table(string(analysisLabel), string(labelA), string(labelB), ...
        metricNames(j), valuesA(j), valuesB(j), valuesA(j)-valuesB(j), ...
        ci(1), ci(2), min(1, pBoot), NaN, ...
        'VariableNames', {'analysis','modelA','modelB','metric','estimateA', ...
        'estimateB','difference','ciLower','ciUpper','pairedBootstrapP','delongP'});
end

[~, pDeLong] = delongTest(y, scoreA, scoreB);
T.delongP(T.metric == "auc") = pDeLong;
end
