function figureTransitionUncertainty(P)
%FIGURETRANSITIONUNCERTAINTY Per-source mean and stochastic-fit variability.

T = readtable(fullfile(P.sourceDir, 'transition_variability_by_source_sample.csv'), ...
    'TextType', 'string', 'VariableNamingRule', 'preserve');
transitions = ["M11","M14","M31"];
abundanceVars = ["sourceAbundanceX1","sourceAbundanceX1","sourceAbundanceX3"];
meanVars = ["meanM11","meanM14","meanM31"];
sdVars = ["sdM11","sdM14","sdM31"];

fig = figure('Position', [60 60 1500 850]);
tl = tiledlayout(fig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
for j = 1:3
    x = T.(abundanceVars(j));
    mu = T.(meanVars(j));
    sd = T.(sdVars(j));
    valid = isfinite(x) & isfinite(mu) & isfinite(sd);
    x = x(valid); mu = mu(valid); sd = sd(valid);
    low = x <= 0.01;

    ax = nexttile(tl, j); hold(ax, 'on');
    shadeLowRegion(ax, 0.01);
    errorbar(ax, x(~low), mu(~low), sd(~low), 'o', ...
        'Color', [0.25 0.50 0.75 0.35], 'MarkerFaceColor', [0.25 0.50 0.75], ...
        'MarkerEdgeColor', 'none', 'MarkerSize', 4, 'CapSize', 0);
    errorbar(ax, x(low), mu(low), sd(low), 'o', ...
        'Color', [0.80 0.25 0.20 0.55], 'MarkerFaceColor', [0.80 0.25 0.20], ...
        'MarkerEdgeColor', 'none', 'MarkerSize', 4, 'CapSize', 0);
    hold(ax, 'off'); xlim(ax, [0 1]); grid(ax, 'on'); box(ax, 'off');
    xlabel(ax, sourceLabel(transitions(j)));
    ylabel(ax, sprintf('%s mean \x00B1 SD across fits', transitions(j)));
    title(ax, sprintf('%s estimates (n=%d)', transitions(j), numel(x)));

    ax = nexttile(tl, j+3); hold(ax, 'on');
    shadeLowRegion(ax, 0.01);
    scatter(ax, x(~low), sd(~low), 34, [0.25 0.50 0.75], 'filled', ...
        'MarkerFaceAlpha', 0.65);
    scatter(ax, x(low), sd(low), 40, [0.80 0.25 0.20], 'filled', ...
        'MarkerFaceAlpha', 0.8);
    [rho, p] = corr(x, sd, 'Type', 'Spearman');
    hold(ax, 'off'); xlim(ax, [0 1]); grid(ax, 'on'); box(ax, 'off');
    xlabel(ax, sourceLabel(transitions(j)));
    ylabel(ax, sprintf('%s SD across fits', transitions(j)));
    title(ax, sprintf('Spearman \\rho=%.2f, p=%.3g', rho, p));
end
title(tl, ['Stochastic-fit uncertainty of reported transitions ' ...
    '(red: source-state abundance \leq 1%)'], 'FontSize', 15, 'FontWeight', 'bold');
exportResponseFigure(fig, P, 'Figure_R4_transition_identifiability');
close(fig);
end

function shadeLowRegion(ax, threshold)
xline(ax, threshold, ':', '1%', 'Color', [0.75 0.20 0.20], ...
    'LabelVerticalAlignment', 'bottom', 'HandleVisibility', 'off');
end

function s = sourceLabel(transition)
if transition == "M31"
    s = 'Observed source-state abundance (x_3)';
else
    s = 'Observed source-state abundance (x_1)';
end
end
