function make_deltaBIC_selected_auc_ba_comparison_figure()
%MAKE_DELTABIC_SELECTED_AUC_BA_COMPARISON_FIGURE
% One figure for BCR::ABL1 and one for MRD. Each has Ridge, Lasso, and
% Elastic net. The winning subset is refit on every outer-training fold.
% Intervals are bootstrap 95% CIs of that subset's held-out scores.

thisDir = fileparts(mfilename('fullpath'));
figureRoot = fileparts(thisDir);
addpath(figureRoot);

markovFile = findSummary(figureRoot, 'markov__deltaBIC2__final_subset_summary.csv');
flowFile = findSummary(figureRoot, 'flow__deltaBIC2__final_subset_summary.csv');

M = readtable(markovFile, 'TextType', 'string', 'VariableNamingRule', 'preserve');
F = readtable(flowFile, 'TextType', 'string', 'VariableNamingRule', 'preserve');

models = {
    "ridge", "Ridge"
    "lasso", "Lasso"
    "elasticnet", "Elastic net"
    };
outcomes = {
    "BCR", "BCR::ABL1", "Figure_deltaBIC_BCR_outer_auc_ba"
    "MRDflow", "MRD", "Figure_deltaBIC_MRD_outer_auc_ba"
    };
viewOrder = ["flow", "markov"];
viewLabels = ["Flow proportion", "Markov"];
viewColors = [
    65 166  89
    51 115 191
] / 255;
metricNames = ["Balanced accuracy", "AUC"];

figureDir = fullfile(figureRoot, 'outputs', 'response_figures');
if ~exist(figureDir, 'dir')
    mkdir(figureDir);
end

picked = table();
for e = 1:size(outcomes, 1)
    outcome = outcomes{e, 1};
    outcomeLabel = outcomes{e, 2};
    stemBase = char(outcomes{e, 3});

for labeled = [1 0]
    if labeled
        stem = stemBase;
    else
        stem = [stemBase '_notext'];
    end

    if labeled
        figWidth = 2100;
        figHeight = 1180;
        metricGap = 1.45;
        rowGap = 3.45;
        yLim = [0.15 7.20];
        xLim = [0.35 1.00];
    else
        figWidth = 2100;
        figHeight = 760;
        metricGap = 0.85;
        rowGap = 2.35;
        yLim = [0.35 5.15];
        xLim = [0.35 1.00];
    end
    fig = figure('Color', 'w', 'Units', 'pixels', ...
        'Position', [40 40 figWidth figHeight], 'Visible', 'off', ...
        'MenuBar', 'none', 'ToolBar', 'none');
    tl = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

    for j = 1:size(models, 1)
        model = models{j, 1};
        panelTitle = models{j, 2};

        ax = nexttile(tl);
        hold(ax, 'on');
        yTicks = zeros(1, numel(viewOrder));
        for v = 1:numel(viewOrder)
            if viewOrder(v) == "markov"
                S = M;
            else
                S = F;
            end
            row = S.outcome == outcome & S.model == model;
            assert(nnz(row) == 1, 'Expected one row for %s %s %s', ...
                viewOrder(v), outcome, model);
            R = S(row, :);
            color = viewColors(v, :);

            subsetName = R.recommendedFinalSubset;
            subsetCount = R.outerSelectionsOfRecommendedSubset;
            yBase = 1 + (v - 1) * rowGap;
            yTicks(v) = yBase + 0.55 * metricGap;
            estimates = [R.nestedBalancedAccuracy, R.nestedAUC];
            lowers = [R.balancedAccuracyCILower, R.aucCILower];
            uppers = [R.balancedAccuracyCIUpper, R.aucCIUpper];
            for t = 1:2
                y = yBase + (t - 1) * metricGap;
                drawInterval(ax, estimates(t), lowers(t), uppers(t), y, color);
                if labeled
                    text(ax, estimates(t), y + 0.32, ...
                        sprintf('%.2f  [%.2f, %.2f]', ...
                        estimates(t), lowers(t), uppers(t)), ...
                        'HorizontalAlignment', 'center', 'FontName', 'Arial', ...
                        'FontSize', 22, 'FontWeight', 'bold', 'Color', 'k');
                    text(ax, estimates(t), y + 0.82, metricNames(t), ...
                        'HorizontalAlignment', 'center', 'FontName', 'Arial', ...
                        'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k');
                end
            end

            if ~labeled
                continue
            end
            picked = [picked; table( ...
                repmat(viewOrder(v), 2, 1), ...
                repmat(string(outcome), 2, 1), ...
                repmat(string(model), 2, 1), ...
                repmat(subsetName, 2, 1), ...
                repmat(subsetCount, 2, 1), ...
                repmat(R.numberOuterSelections, 2, 1), ...
                metricNames(:), estimates(:), lowers(:), uppers(:), ...
                'VariableNames', {'view','outcome','model','selectedFeatures', ...
                'outerSelectionsOfSubset','nOuterFolds','metric', ...
                'estimate','ciLower','ciUpper'})]; %#ok<AGROW>
        end

        ax.XLim = xLim;
        ax.YLim = yLim;
        ax.YDir = 'reverse';
        ax.YTick = yTicks;
        ax.YTickLabel = viewLabels;
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
        if labeled
            ax.XLabel.FontSize = 16;
            ax.XLabel.FontWeight = 'bold';
            title(ax, panelTitle, 'FontWeight', 'bold', 'FontSize', 24, 'Color', 'k');
        else
            ax.XLabel.FontSize = 13;
            ax.XLabel.FontWeight = 'normal';
        end
    end

    if labeled
        title(tl, [outcomeLabel '  |  outer-test AUC and balanced accuracy'], ...
            'FontName', 'Arial', 'FontSize', 26, 'FontWeight', 'bold', ...
            'Interpreter', 'none');
    end

    axs = findall(fig, 'Type', 'axes');
    for i = 1:numel(axs)
        if isprop(axs(i), 'Toolbar') && ~isempty(axs(i).Toolbar)
            axs(i).Toolbar.Visible = 'off';
        end
        disableDefaultInteractivity(axs(i));
    end
    exportgraphics(fig, fullfile(figureDir, [stem '.png']), 'Resolution', 300);
    exportgraphics(fig, fullfile(figureDir, [stem '.pdf']), 'ContentType', 'vector');
    close(fig);
    fprintf('Wrote %s\n', fullfile(figureDir, [stem '.png']));
end
end

writetable(picked, fullfile(figureDir, 'Figure_deltaBIC_BCR_MRD_outer_auc_ba.csv'));
end

function summaryFile = findSummary(figureRoot, name)
candidates = {
    fullfile(figureRoot, 'outputs_10fold', 'delta_bic_selection', name)
    fullfile(figureRoot, 'outputs', 'delta_bic_selection', name)
    };
for i = 1:numel(candidates)
    if exist(candidates{i}, 'file') == 2
        summaryFile = candidates{i};
        return
    end
end
error('Missing %s. Run the nested search and select_final_subset_nested_delta_bic first.', name);
end

function drawInterval(ax, point, lower, upper, y, color)
line(ax, [lower upper], [y y], 'Color', color, 'LineWidth', 2.2);
line(ax, [lower lower], [y - 0.08 y + 0.08], 'Color', color, 'LineWidth', 1.6);
line(ax, [upper upper], [y - 0.08 y + 0.08], 'Color', color, 'LineWidth', 1.6);
plot(ax, point, y, 'd', 'MarkerSize', 9, ...
    'MarkerFaceColor', color, 'MarkerEdgeColor', 'w', 'LineWidth', 1.0);
end
