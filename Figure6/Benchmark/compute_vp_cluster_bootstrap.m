function files = compute_vp_cluster_bootstrap(scoreDir)
%COMPUTE_VP_CLUSTER_BOOTSTRAP Patient-cluster bootstrap of VP mean-OOF scores.
%
% A source block is the 100 virtual-patient rows generated for one patient
% sample (VirtualPatient restarts at 1). BCR has 101 blocks and MRD has 104.
% Some source ID values appear in more than one block, so the resampled unit
% is the block, not the raw ID string. A resample draws those blocks with
% replacement, within class, and keeps every VP row of each drawn block.
% The point estimate remains the metric on all VP rows. Intervals are the
% 2.5th and 97.5th percentiles of 2,000 resamples.

thisDir = fileparts(mfilename('fullpath'));
figureRoot = fileparts(thisDir);
if nargin < 1 || isempty(scoreDir)
    scoreDir = fullfile(figureRoot, 'outputs', 'specimen_folds');
end
completeVPFile = fullfile(figureRoot, 'Data', 'PB_BM_merged_with_MRDflow_BCR.csv');
assert(isfile(completeVPFile), 'Missing VP table: %s', completeVPFile);
addpath(figureRoot);
addpath(fullfile(figureRoot, 'Functions'));

nBoot = 2000;
seed = 1729 + 1101;
cohorts = ["BCR","MRDflow"];
files = struct();

T = readtable(completeVPFile, 'VariableNamingRule', 'preserve', 'TextType', 'string');
id = string(T.ID);
vp = double(T.VirtualPatient);
sourceBlock = cumsum(vp == 1);
mrd = strtrim(string(T{:, 'MRD(flow)'}));
bcr = strtrim(string(T{:, 'BCR::ABL1/like'}));

for i = 1:numel(cohorts)
    if cohorts(i) == "BCR"
        keep = ismember(bcr, ["Negative","BCR::ABL1"]);
        y = double(bcr(keep) == "BCR::ABL1");
        expectedClusters = 101;
    else
        keep = ismember(mrd, ["Negative","Positive (< 10^-3)","Positive (> 10^-3)"]);
        y = double(mrd(keep) == "Positive (> 10^-3)");
        expectedClusters = 104;
    end
    blockID = sourceBlock(keep);
    sourceID = id(keep);
    [patientBlock, ~, cluster] = unique(blockID, 'stable');
    nClusters = numel(patientBlock);
    assert(nClusters == expectedClusters, ...
        '%s expected %d source blocks, found %d.', ...
        cohorts(i), expectedClusters, nClusters);
    rowsByCluster = cell(nClusters, 1);
    clusterClass = zeros(nClusters, 1);
    for c = 1:nClusters
        rows = find(cluster == c);
        rowsByCluster{c} = rows;
        classes = unique(y(rows));
        assert(isscalar(classes), 'Source block %d has mixed outcomes.', patientBlock(c));
        assert(numel(unique(sourceID(rows))) == 1, ...
            'Source block %d contains more than one ID.', patientBlock(c));
        clusterClass(c) = classes;
        assert(numel(rows) == 100, ...
            'Source block %d has %d VP rows; expected 100.', patientBlock(c), numel(rows));
    end
    posClusters = find(clusterClass == 1);
    negClusters = find(clusterClass == 0);
    assert(numel(y) == 100 * nClusters, 'VP row count does not match 100 per patient.');

    scores = loadCohortScores(scoreDir, cohorts(i), y);
    modelNames = string(scores.Properties.VariableNames);
    summary = clusterBootstrap(y, scores{:, :}, modelNames, ...
        rowsByCluster, posClusters, negClusters, nBoot, seed, ...
        "VP__" + cohorts(i) + "__markov");
    outFile = fullfile(scoreDir, ...
        char("vp_" + lower(cohorts(i)) + "__markov_cluster_bootstrap_95CI.csv"));
    writetable(summary, outFile);
    files.(char(cohorts(i))) = string(outFile);
    fprintf('Wrote %s\n', outFile);
end
end

function scores = loadCohortScores(scoreDir, cohort, y)
baseFile = fullfile(scoreDir, ...
    char("vp_" + lower(cohort) + "__markov_mean_OOF_scores.csv"));
ridgeFile = fullfile(scoreDir, ...
    char("vp_" + lower(cohort) + "__markov_ridge_mean_OOF_scores.csv"));
base = readtable(baseFile);
ridge = readtable(ridgeFile);
assert(isequal(base.outcome, y) && isequal(ridge.outcome, y), ...
    'Mean-OOF outcome column does not match the VP cohort filter.');
assert(isequal(base.row, ridge.row), 'Ridge and Lasso mean-OOF row order differ.');
scores = table(base.lasso_meanOOF, base.elasticnet_meanOOF, ridge.ridge_meanOOF, ...
    'VariableNames', {'lasso','elasticnet','ridge'});
end

function T = clusterBootstrap(y, scores, modelNames, rowsByCluster, ...
        posClusters, negClusters, nBoot, seed, analysisLabel)
y = double(y(:));
scores = double(scores);
modelNames = string(modelNames(:));
nModels = size(scores, 2);
nPos = numel(posClusters);
nNeg = numel(negClusters);
n = nPos + nNeg;
blockRows = 100;
assert(all(cellfun(@numel, rowsByCluster) == blockRows), ...
    'Cluster bootstrap expects 100 VP rows per patient.');

rng(seed, 'twister');
metricNames = ["accuracy","balancedAccuracy","sensitivity","specificity","auc"];
boot = nan(nBoot, numel(metricNames), nModels);
for b = 1:nBoot
    draw = [posClusters(randi(nPos, nPos, 1)); ...
            negClusters(randi(nNeg, nNeg, 1))];
    idx = zeros(n * blockRows, 1);
    for k = 1:n
        idx((k - 1) * blockRows + (1:blockRows)) = rowsByCluster{draw(k)};
    end
    yb = y(idx);
    for m = 1:nModels
        M = metrics(yb, scores(idx, m));
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
            n, nPos, nNeg, "patient_cluster", ...
            'VariableNames', {'analysis','model','metric','estimate','ciLower','ciUpper', ...
            'n','nPositive','nNegative','nPatients','nPositivePatients', ...
            'nNegativePatients','bootstrap'});
    end
end
end
