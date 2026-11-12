function saveNestedResult(R, cfg, stem)
%SAVENESTEDRESULT Save repeat metrics, CIs, stability, OOF scores, and MAT.

stem = matlab.lang.makeValidName(char(stem));
writetable(R.repeatMetrics, fullfile(cfg.outputDir, [stem '_repeat_metrics.csv']));
writetable(R.selectionTable, fullfile(cfg.outputDir, [stem '_selection_stability.csv']));

summary = bootstrapSummary(R.y, R.meanOOFScores, R.modelNames, ...
    cfg.bootstrapIterations, cfg.randomSeed + 901, R.analysisLabel);
writetable(summary, fullfile(cfg.outputDir, [stem '_bootstrap_95CI.csv']));

oof = table((1:numel(R.y))', R.y, 'VariableNames', {'row','outcome'});
for m = 1:numel(R.modelNames)
    oof.(matlab.lang.makeValidName(R.modelNames(m) + "_meanOOF")) = R.meanOOFScores(:, m);
end
writetable(oof, fullfile(cfg.outputDir, [stem '_mean_OOF_scores.csv']));
save(fullfile(cfg.outputDir, [stem '.mat']), 'R', 'summary', 'cfg', '-v7.3');
end
