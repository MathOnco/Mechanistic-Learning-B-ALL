function make_balanced_mrd_outer_auc_ba_figure()
%MAKE_BALANCED_MRD_OUTER_AUC_BA_FIGURE Winner outer-test performance.
%
% Each point uses the Delta-BIC winner chosen inside that outer fold.
% Markov values are the published original-data nested results. Balanced VP
% values use the same winner rule on the one balanced flow-MRD draw. The VP
% interval resamples source specimens, keeping every selected virtual-patient
% row from a drawn specimen (2,000 resamples).

thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, 'Functions'));
detailDir = fullfile(thisDir, 'outputs');
figureDir = fullfile(detailDir, 'figures');
if ~isfolder(figureDir), mkdir(figureDir); end

markovFile = fullfile(thisDir, 'Data', 'markov__deltaBIC2__final_subset_summary.csv');
Markov = readtable(markovFile, 'TextType', 'string', 'VariableNamingRule', 'preserve');
Markov = Markov(Markov.outcome == "MRDflow", :);

models = ["ridge","lasso","elasticnet"];
modelLabels = ["Ridge","Lasso","Elastic net"];
representations = ["Markov","Balanced VP"];
colors = [51 115 191; 0 140 130] / 255;
metricNames = ["Balanced accuracy","AUC"];
nBoot = 2000;

VP = table();
for m = 1:numel(models)
    stem = "mrdflow__markov__" + models(m);
    saved = load(fullfile(detailDir, char(stem + "__result.mat")), 'R', 'C');
    selected = readtable(fullfile(figureDir, ...
        char(stem + "__deltaBIC2__fold_selections.csv")), ...
        'TextType', 'string', 'VariableNamingRule', 'preserve');
    score = selectedOuterScore(saved.R, selected);
    point = metrics(saved.C.y, score);
    [lower, upper] = specimenBootstrap(saved.C.y, score, saved.C.group, ...
        nBoot, 1729 + m);
    VP = [VP; table("markov", "MRDflow", models(m), point.n, ...
        point.balancedAccuracy, lower(1), upper(1), ...
        point.auc, lower(2), upper(2), ...
        'VariableNames', {'view','outcome','model','n', ...
        'nestedBalancedAccuracy','balancedAccuracyCILower','balancedAccuracyCIUpper', ...
        'nestedAUC','aucCILower','aucCIUpper'})]; %#ok<AGROW>
end
writetable(VP, fullfile(figureDir, 'balanced_mrd_vp_deltaBIC2_outer_auc_ba.csv'));

fig = figure('Color', 'w', 'Units', 'pixels', ...
    'Position', [40 40 2100 1180], 'Visible', 'off', ...
    'MenuBar', 'none', 'ToolBar', 'none');
tl = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
metricGap = 1.45;
rowGap = 3.45;

for j = 1:numel(models)
    ax = nexttile(tl);
    hold(ax, 'on');
    yTicks = zeros(1, numel(representations));
    for v = 1:numel(representations)
        if representations(v) == "Markov"
            row = Markov.model == models(j);
            R = Markov(row, :);
        else
            row = VP.model == models(j);
            R = VP(row, :);
        end
        assert(height(R) == 1, 'Expected one %s %s row.', representations(v), models(j));
        yBase = 1 + (v - 1) * rowGap;
        yTicks(v) = yBase + 0.55 * metricGap;
        estimates = [R.nestedBalancedAccuracy, R.nestedAUC];
        lowers = [R.balancedAccuracyCILower, R.aucCILower];
        uppers = [R.balancedAccuracyCIUpper, R.aucCIUpper];
        for t = 1:2
            y = yBase + (t - 1) * metricGap;
            drawInterval(ax, estimates(t), lowers(t), uppers(t), y, colors(v, :));
            text(ax, estimates(t), y + 0.32, ...
                sprintf('%.2f  [%.2f, %.2f]', estimates(t), lowers(t), uppers(t)), ...
                'HorizontalAlignment', 'center', 'FontName', 'Arial', ...
                'FontSize', 22, 'FontWeight', 'bold', 'Color', 'k');
            text(ax, estimates(t), y + 0.82, metricNames(t), ...
                'HorizontalAlignment', 'center', 'FontName', 'Arial', ...
                'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k');
        end
    end
    ax.XLim = [0.35 1.00];
    ax.YLim = [0.15 7.20];
    ax.YDir = 'reverse';
    ax.YTick = yTicks;
    ax.YTickLabel = representations;
    ax.XTick = 0.4:0.1:1.0;
    ax.XGrid = 'on';
    ax.YGrid = 'off';
    ax.GridAlpha = 0.16;
    ax.TickLength = [0 0];
    ax.FontName = 'Arial';
    ax.Box = 'on';
    xlabel(ax, {'Nested-CV outer-test performance', ...
        'Bootstrap 95% CI, 2,000 resamples'});
    ax.FontSize = 22;
    ax.FontWeight = 'bold';
    ax.XAxis.FontSize = 22;
    ax.XAxis.FontWeight = 'bold';
    ax.YAxis.FontSize = 18;
    ax.YAxis.FontWeight = 'bold';
    ax.XLabel.FontSize = 16;
    ax.XLabel.FontWeight = 'bold';
    title(ax, char(modelLabels(j)), 'FontWeight', 'bold', 'FontSize', 24, 'Color', 'k');
end

title(tl, 'MRD  |  outer-test AUC and balanced accuracy', ...
    'FontName', 'Arial', 'FontSize', 26, 'FontWeight', 'bold', ...
    'Interpreter', 'none');
axs = findall(fig, 'Type', 'axes');
for i = 1:numel(axs)
    if isprop(axs(i), 'Toolbar') && ~isempty(axs(i).Toolbar)
        axs(i).Toolbar.Visible = 'off';
    end
    disableDefaultInteractivity(axs(i));
end

stemOut = 'Figure_balanced_MRD_VP_vs_Markov_outer_auc_ba';
exportgraphics(fig, fullfile(figureDir, [stemOut '.png']), 'Resolution', 300);
close(fig);
fprintf('Wrote %s\n', fullfile(figureDir, [stemOut '.png']));
end

function score = selectedOuterScore(R, selected)
n = size(R.oofScores, 1);
score = nan(n, 1);
for i = 1:height(selected)
    r = selected.repeat(i);
    f = selected.outerFold(i);
    k = selected.k(i);
    testRows = R.outerFoldIndex(:, r) == f;
    score(testRows) = R.oofScores(testRows, r, k);
end
assert(all(isfinite(score)), 'A row is missing its winner outer-test score.');
end

function [lower, upper] = specimenBootstrap(y, score, group, iterations, seed)
y = double(y(:));
score = double(score(:));
group = double(group(:));
[groups, ~, groupIndex] = unique(group, 'stable');
rowsByGroup = cell(numel(groups), 1);
classOfGroup = zeros(numel(groups), 1);
for g = 1:numel(groups)
    rows = find(groupIndex == g);
    rowsByGroup{g} = rows;
    labels = unique(y(rows));
    assert(isscalar(labels), 'Specimen %g has mixed outcomes.', groups(g));
    classOfGroup(g) = labels;
end
pos = find(classOfGroup == 1);
neg = find(classOfGroup == 0);
rng(seed, 'twister');
boot = nan(iterations, 2);
for b = 1:iterations
    draw = [pos(randi(numel(pos), numel(pos), 1)); ...
            neg(randi(numel(neg), numel(neg), 1))];
    idx = vertcat(rowsByGroup{draw});
    M = metrics(y(idx), score(idx));
    boot(b, :) = [M.balancedAccuracy, M.auc];
end
lower = prctile(boot, 2.5, 1);
upper = prctile(boot, 97.5, 1);
end

function drawInterval(ax, point, lower, upper, y, color)
line(ax, [lower upper], [y y], 'Color', color, 'LineWidth', 2.2);
line(ax, [lower lower], [y - 0.08 y + 0.08], 'Color', color, 'LineWidth', 1.6);
line(ax, [upper upper], [y - 0.08 y + 0.08], 'Color', color, 'LineWidth', 1.6);
plot(ax, point, y, 'd', 'MarkerSize', 9, ...
    'MarkerFaceColor', color, 'MarkerEdgeColor', 'w', 'LineWidth', 1.0);
end
