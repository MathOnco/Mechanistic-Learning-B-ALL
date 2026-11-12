function analysis_01_nested_classifier_benchmark()
%ANALYSIS_01_NESTED_CLASSIFIER_BENCHMARK Nested classifier benchmarking.
% Compares ridge, lasso, elastic net, and random forest under nested CV
% with lambda tuning, bootstrap intervals, and pairwise model comparisons.

clc; close all;
thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, 'Functions'));
cfg = benchmark_config();

outcomes = ["MRDflow", "BCR"];
allSummary = table();
allComparisons = table();

for o = 1:numel(outcomes)
    outcome = outcomes(o);
    C = loadCohort(cfg, outcome);
    fprintf('\n%s: N=%d, positive=%d, negative=%d\n', ...
        outcome, numel(C.y), sum(C.y == 1), sum(C.y == 0));

    results = cell(numel(cfg.views), 1);
    for v = 1:numel(cfg.views)
        viewName = cfg.views(v);
        [X, featureNames] = encodeView(C, viewName);
        label = outcome + "__" + viewName;
        fprintf('  Nested CV: %s (%d predictors)\n', label, size(X, 2));
        R = nestedCV(X, C.y, featureNames, cfg.models, cfg, label);
        results{v} = R;
        stem = lower(char(label));
        saveNestedResult(R, cfg, stem);

        S = bootstrapSummary(C.y, R.meanOOFScores, R.modelNames, ...
            cfg.bootstrapIterations, cfg.randomSeed + 1000*o + v, label);
        allSummary = [allSummary; S]; %#ok<AGROW>
    end

    markovIdx = find(cfg.views == "markov", 1);
    referenceViews = ["flow", "clinical_genomic"];
    for rv = 1:numel(referenceViews)
        refIdx = find(cfg.views == referenceViews(rv), 1);
        if isempty(markovIdx) || isempty(refIdx), continue; end
        A = results{markovIdx};
        B = results{refIdx};
        for m = 1:numel(cfg.models)
            comparisonLabel = outcome + "__" + cfg.models(m) + ...
                "__markov_vs_" + referenceViews(rv);
            Tcmp = pairedModelComparison(C.y, ...
                A.meanOOFScores(:, m), B.meanOOFScores(:, m), ...
                "markov_" + cfg.models(m), referenceViews(rv) + "_" + cfg.models(m), ...
                cfg.bootstrapIterations, cfg.randomSeed + 5000*o + 100*rv + m, ...
                comparisonLabel);
            allComparisons = [allComparisons; Tcmp]; %#ok<AGROW>
        end
    end
end

writetable(allSummary, fullfile(cfg.outputDir, 'benchmark_all_bootstrap_95CI.csv'));
writetable(allComparisons, fullfile(cfg.outputDir, 'benchmark_paired_comparisons.csv'));

% Publication-ready compact AUC plot.
aucRows = allSummary(allSummary.metric == "auc", :);
labels = aucRows.analysis + " | " + aucRows.model;
figure('Color', 'w', 'Position', [100 100 1100 max(500, 22*height(aucRows))]);
y = (1:height(aucRows))';
errorbar(aucRows.estimate, y, aucRows.estimate-aucRows.ciLower, ...
    aucRows.ciUpper-aucRows.estimate, 'horizontal', 'o', ...
    'LineWidth', 1.2, 'MarkerFaceColor', [0.2 0.45 0.75]);
xlim([0 1]); ylim([0 height(aucRows)+1]); grid on; box off;
yticks(y); yticklabels(labels); set(gca, 'YDir', 'reverse', 'TickLabelInterpreter', 'none');
xlabel('Mean repeated-OOF AUC (stratified-bootstrap 95% CI)');
exportgraphics(gcf, fullfile(cfg.outputDir, 'benchmark_AUC_95CI.png'), 'Resolution', 300);

fprintf('\nSaved nested benchmark outputs to %s\n', cfg.outputDir);
end
