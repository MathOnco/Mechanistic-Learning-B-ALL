function figureVPAugmentation(P)
%FIGUREVPAUGMENTATION Performance versus generated observations per source.

summaryFile = fullfile(P.sourceDir, ...
    'VP_augmentation_regime_summary_specimen_folds.csv');
if ~isfile(summaryFile)
    summaryFile = fullfile(P.sourceDir, 'VP_augmentation_regime_summary.csv');
end
T = readtable(summaryFile, 'TextType', 'string', 'VariableNamingRule', 'preserve');
outcomes = ["MRDflow","BCR"];
fig = figure('Position', [100 100 1150 480]);
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
models = unique(T.model, 'stable');
colors = lines(numel(models));
for o = 1:2
    ax = nexttile(tl); hold(ax, 'on');
    for m = 1:numel(models)
        S = sortrows(T(T.outcome == outcomes(o) & T.model == models(m), :), ...
            'virtualPatientsPerSource');
        errorbar(ax, S.virtualPatientsPerSource, S.meanAUC, ...
            S.meanAUC-S.aucLower, S.aucUpper-S.meanAUC, '-o', ...
            'LineWidth', 1.8, 'Color', colors(m,:), 'MarkerFaceColor', colors(m,:), ...
            'DisplayName', prettyModel(models(m)));
    end
    hold(ax, 'off'); ylim(ax, [0 1]); grid(ax, 'on'); box(ax, 'off');
    xticks(ax, unique(T.virtualPatientsPerSource));
    xlabel(ax, 'Generated observations sampled per source sample');
    ylabel(ax, 'Nested-CV AUC across VP resamples');
    legend(ax, 'Location', 'best');
    if outcomes(o) == "MRDflow", title(ax, 'Flow MRD'); else, title(ax, 'BCR::ABL1'); end
end
title(tl, ['Outer-test performance and sensitivity analysis using VP ' ...
    'to predict BCR::ABL1 and MRD'], ...
    'FontSize', 14, 'FontWeight', 'bold');
exportResponseFigure(fig, P, 'Figure_R8_VP_augmentation_regime');
close(fig);
end

function s = prettyModel(s)
s = replace(s, "elasticnet", "Elastic net");
s = replace(s, "randomforest", "Random forest");
s = replace(s, "ridge", "Ridge");
s = replace(s, "lasso", "Lasso");
end
