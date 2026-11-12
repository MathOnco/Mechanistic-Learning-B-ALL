function outputs = select_final_subset_nested_delta_bic( ...
    sourceDir, view, deltaBICThreshold, bootstrapIterations)
%SELECT_FINAL_SUBSET_NESTED_DELTA_BIC Foldwise BIC-screened subset selection.
%
% outputs = select_final_subset_nested_delta_bic()
% outputs = select_final_subset_nested_delta_bic(sourceDir, view, 2, 2000)
%
% This is a post-processing analysis. It does not rerun the exhaustive
% searches. Within each outer fold it:
%   1. finds the minimum training BIC among the best-k candidates;
%   2. retains candidates having Delta BIC <= deltaBICThreshold;
%   3. chooses the largest INNER-CV selection score among those candidates;
%   4. breaks exact score ties toward smaller k, then smaller BIC/mask;
%   5. uses that row's already-saved outer-test predictions.
%
% The modal exact subset across outer folds is reported as the recommended
% descriptive full-cohort subset. Generalization performance comes from the
% fold-specific selected models, not from refitting the modal subset.

thisDir = fileparts(mfilename('fullpath'));
addpath(thisDir);
if nargin < 1 || isempty(sourceDir)
    sourceDir = fullfile(thisDir, 'outputs');
end
if nargin < 2 || isempty(view), view = "markov"; end
if nargin < 3 || isempty(deltaBICThreshold), deltaBICThreshold = 2; end
if nargin < 4 || isempty(bootstrapIterations), bootstrapIterations = 2000; end

sourceDir = char(string(sourceDir));
view = lower(string(view));
validateattributes(deltaBICThreshold, {'numeric'}, ...
    {'scalar','real','finite','nonnegative'});
validateattributes(bootstrapIterations, {'numeric'}, ...
    {'scalar','integer','positive'});

pattern = char("*__" + view + "__*__outer_fold_details.csv");
files = dir(fullfile(sourceDir, pattern));
if isempty(files)
    error('deltaBIC:NoFiles', ...
        'No outer-fold files matching %s were found in %s.', pattern, sourceDir);
end

tag = deltaTag(deltaBICThreshold);
outputDir = fullfile(sourceDir, 'delta_bic_selection');
if ~isfolder(outputDir), mkdir(outputDir); end

allSummary = table();
allPerformance = table();
allFeatureFrequency = table();

for i = 1:numel(files)
    detailFile = fullfile(files(i).folder, files(i).name);
    stem = erase(string(files(i).name), "__outer_fold_details.csv");
    resultFile = fullfile(sourceDir, char(stem + "__result.mat"));
    if ~isfile(resultFile)
        error('deltaBIC:MissingResult', ...
            'The result MAT file required for patient-level OOF scores is missing: %s', ...
            resultFile);
    end

    D = readtable(detailFile, 'TextType', 'string', ...
        'VariableNamingRule', 'preserve');
    validateDetailTable(D, detailFile);
    selected = chooseWithinOuterFolds(D, deltaBICThreshold);

    M = load(resultFile, 'R', 'C');
    [oof, metricRow] = assembleSelectedOOF(M.R, M.C, selected, ...
        bootstrapIterations, 1729 + i);
    [summaryRow, featureFrequency] = summarizeSelections( ...
        selected, metricRow, M.C.featureNames, deltaBICThreshold);

    writetable(selected, fullfile(outputDir, ...
        char(stem + "__" + tag + "__fold_selections.csv")));
    writetable(oof, fullfile(outputDir, ...
        char(stem + "__" + tag + "__selected_oof_scores.csv")));
    writetable(featureFrequency, fullfile(outputDir, ...
        char(stem + "__" + tag + "__feature_frequency.csv")));

    allSummary = [allSummary; summaryRow]; %#ok<AGROW>
    allPerformance = [allPerformance; metricRow]; %#ok<AGROW>
    allFeatureFrequency = [allFeatureFrequency; featureFrequency]; %#ok<AGROW>
end

allSummary = sortrows(allSummary, {'outcome','model'});
allPerformance = sortrows(allPerformance, {'outcome','model'});
allFeatureFrequency = sortrows(allFeatureFrequency, ...
    {'outcome','model','selectionFrequency'}, {'ascend','ascend','descend'});

summaryFile = fullfile(outputDir, ...
    char(view + "__" + tag + "__final_subset_summary.csv"));
performanceFile = fullfile(outputDir, ...
    char(view + "__" + tag + "__performance_95CI.csv"));
featureFile = fullfile(outputDir, ...
    char(view + "__" + tag + "__all_feature_frequencies.csv"));
writetable(allSummary, summaryFile);
writetable(allPerformance, performanceFile);
writetable(allFeatureFrequency, featureFile);

outputs = struct('summary', string(summaryFile), ...
    'performance', string(performanceFile), ...
    'featureFrequency', string(featureFile), ...
    'outputDir', string(outputDir));
fprintf('Delta-BIC nested selection outputs written to %s\n', outputDir);
end

function validateDetailTable(D, filename)
required = ["view","outcome","model","fixedLambda","repeat", ...
    "outerFold","k","selectedMask","selectedFeatures", ...
    "innerSelectionValue","trainingBIC"];
missing = required(~ismember(required, string(D.Properties.VariableNames)));
if ~isempty(missing)
    error('deltaBIC:MissingColumns', '%s is missing: %s', ...
        filename, strjoin(missing, ', '));
end
if any(~isfinite(D.trainingBIC)) || any(~isfinite(D.innerSelectionValue))
    error('deltaBIC:NonfiniteSelectionValues', ...
        '%s contains nonfinite BIC or inner selection values.', filename);
end
end

function selected = chooseWithinOuterFolds(D, threshold)
pairs = unique(D(:, {'repeat','outerFold'}), 'rows', 'stable');
D.deltaBIC = nan(height(D), 1);
D.minimumTrainingBIC = nan(height(D), 1);
D.numberBICEligible = zeros(height(D), 1);
selected = D([], :);
for i = 1:height(pairs)
    inFold = D.repeat == pairs.repeat(i) & D.outerFold == pairs.outerFold(i);
    F = D(inFold, :);
    minimumBIC = min(F.trainingBIC);
    F.deltaBIC = F.trainingBIC - minimumBIC;
    eligible = F(F.deltaBIC <= threshold + 1e-12, :);
    eligible = sortrows(eligible, ...
        {'innerSelectionValue','k','trainingBIC','selectedMask'}, ...
        {'descend','ascend','ascend','ascend'});
    winner = eligible(1, :);
    winner.minimumTrainingBIC = minimumBIC;
    winner.numberBICEligible = height(eligible);
    selected = [selected; winner]; %#ok<AGROW>
end
selected = sortrows(selected, {'repeat','outerFold'});
end

function [oof, metricRow] = assembleSelectedOOF( ...
    R, C, selected, bootstrapIterations, seed)
n = numel(C.y);
nRepeats = R.meta.nRepeats;
scoresByRepeat = nan(n, nRepeats);
selectedKByRepeat = nan(n, nRepeats);

for i = 1:height(selected)
    r = selected.repeat(i);
    f = selected.outerFold(i);
    k = selected.k(i);
    testRows = R.outerFoldIndex(:, r) == f;
    scoresByRepeat(testRows, r) = R.oofScores(testRows, r, k);
    selectedKByRepeat(testRows, r) = k;
end
if any(~isfinite(scoresByRepeat), 'all')
    error('deltaBIC:IncompleteOOF', ...
        'At least one observation lacks a selected outer-fold prediction.');
end

meanScore = mean(scoresByRepeat, 2);
meanSelectedK = mean(selectedKByRepeat, 2);
oof = table(C.originalIndex, C.patientID, C.sampleID, C.y, ...
    meanSelectedK, meanScore, ...
    'VariableNames', {'originalIndex','patientID','sampleID','outcome', ...
    'meanSelectedK','selectedOOFScore'});
for r = 1:nRepeats
    oof.(sprintf('scoreRepeat%02d', r)) = scoresByRepeat(:, r);
    oof.(sprintf('selectedKRepeat%02d', r)) = selectedKByRepeat(:, r);
end

point = exhaustive.metrics(C.y, meanScore);
[lower, upper] = bootstrapMetrics( ...
    C.y, meanScore, bootstrapIterations, seed);
selectionScore = 0.5 * point.balancedAccuracy + 0.5 * point.auc;
metricRow = table(R.view, R.outcome, R.model, R.lambda, point.n, ...
    point.nPositive, point.nNegative, bootstrapIterations, ...
    point.accuracy, lower(1), upper(1), ...
    point.balancedAccuracy, lower(2), upper(2), ...
    point.auc, lower(3), upper(3), ...
    selectionScore, lower(4), upper(4), ...
    'VariableNames', {'view','outcome','model','fixedLambda','n', ...
    'nPositive','nNegative','bootstrapIterations', ...
    'accuracy','accuracyCILower','accuracyCIUpper', ...
    'balancedAccuracy','balancedAccuracyCILower','balancedAccuracyCIUpper', ...
    'auc','aucCILower','aucCIUpper', ...
    'selectionScore','selectionScoreCILower','selectionScoreCIUpper'});
end

function [lower, upper] = bootstrapMetrics(y, score, iterations, seed)
y = double(y(:));
score = double(score(:));
positive = find(y == 1);
negative = find(y == 0);
if isempty(positive) || isempty(negative)
    error('deltaBIC:OneClass', 'Both outcome classes are required.');
end

rng(seed, 'twister');
boot = nan(iterations, 4);
for b = 1:iterations
    ip = positive(randi(numel(positive), numel(positive), 1));
    in = negative(randi(numel(negative), numel(negative), 1));
    idx = [ip; in];
    M = exhaustive.metrics(y(idx), score(idx));
    combo = 0.5 * M.balancedAccuracy + 0.5 * M.auc;
    boot(b, :) = [M.accuracy, M.balancedAccuracy, M.auc, combo];
end
lower = prctile(boot, 2.5, 1);
upper = prctile(boot, 97.5, 1);
end

function [summary, featureFrequency] = summarizeSelections( ...
    selected, metricRow, featureNames, threshold)
nSelections = height(selected);
masks = unique(selected.selectedMask);
candidates = table();
for i = 1:numel(masks)
    hit = selected.selectedMask == masks(i);
    first = find(hit, 1);
    candidates(i, :) = table(masks(i), selected.k(first), ...
        selected.selectedFeatures(first), sum(hit), ...
        mean(selected.innerSelectionValue(hit)), ...
        mean(selected.trainingBIC(hit)), mean(selected.deltaBIC(hit)), ...
        'VariableNames', {'mask','k','features','count', ...
        'meanInnerSelectionValue','meanTrainingBIC','meanDeltaBIC'});
end
candidates = sortrows(candidates, ...
    {'count','meanInnerSelectionValue','k','meanTrainingBIC','mask'}, ...
    {'descend','descend','ascend','ascend','ascend'});
winner = candidates(1, :);

labels = strings(height(candidates), 1);
for i = 1:height(candidates)
    labels(i) = candidates.features(i) + " (" + ...
        candidates.count(i) + "/" + nSelections + ")";
end
allSubsetFrequencies = strjoin(labels, '; ');

ks = unique(selected.k);
kLabels = strings(numel(ks), 1);
for i = 1:numel(ks)
    count = sum(selected.k == ks(i));
    kLabels(i) = "k=" + ks(i) + " (" + count + "/" + nSelections + ")";
end
kDistribution = strjoin(kLabels, '; ');

selectionRule = "Within each outer fold: Delta BIC <= " + threshold + ...
    "; maximize inner 0.5*balancedAccuracy + 0.5*AUC; " + ...
    "tie -> smaller k";
summary = table(metricRow.view, metricRow.outcome, metricRow.model, ...
    metricRow.fixedLambda, threshold, nSelections, winner.k, winner.mask, ...
    winner.features, winner.count, winner.count / nSelections, ...
    kDistribution, allSubsetFrequencies, selectionRule, ...
    metricRow.accuracy, metricRow.accuracyCILower, metricRow.accuracyCIUpper, ...
    metricRow.balancedAccuracy, metricRow.balancedAccuracyCILower, ...
    metricRow.balancedAccuracyCIUpper, metricRow.auc, ...
    metricRow.aucCILower, metricRow.aucCIUpper, ...
    metricRow.selectionScore, metricRow.selectionScoreCILower, ...
    metricRow.selectionScoreCIUpper, ...
    'VariableNames', {'view','outcome','model','fixedLambda', ...
    'deltaBICThreshold','numberOuterSelections','recommendedFinalK', ...
    'recommendedFinalMask','recommendedFinalSubset', ...
    'outerSelectionsOfRecommendedSubset','recommendedSubsetFrequency', ...
    'selectedKDistribution','allSelectedSubsetFrequencies','selectionRule', ...
    'nestedAccuracy','accuracyCILower','accuracyCIUpper', ...
    'nestedBalancedAccuracy','balancedAccuracyCILower', ...
    'balancedAccuracyCIUpper','nestedAUC','aucCILower','aucCIUpper', ...
    'nestedSelectionScore','selectionScoreCILower', ...
    'selectionScoreCIUpper'});

featureNames = string(featureNames(:));
counts = zeros(numel(featureNames), 1);
for i = 1:height(selected)
    chosen = split(strrep(selected.selectedFeatures(i), " ", ""), ',');
    for j = 1:numel(chosen)
        counts(featureNames == chosen(j)) = counts(featureNames == chosen(j)) + 1;
    end
end
featureFrequency = table(repmat(metricRow.view, numel(featureNames), 1), ...
    repmat(metricRow.outcome, numel(featureNames), 1), ...
    repmat(metricRow.model, numel(featureNames), 1), featureNames, counts, ...
    counts / nSelections, ...
    'VariableNames', {'view','outcome','model','feature','selectionCount', ...
    'selectionFrequency'});
end

function tag = deltaTag(threshold)
tag = "deltaBIC" + replace(string(sprintf('%.3g', threshold)), '.', 'p');
end
