function make_balanced_mrd_fold_overlay_figure()
%MAKE_BALANCED_MRD_FOLD_OVERLAY_FIGURE Inner-CV paths for the balanced MRD VP search.
%
% One panel per model. Each line is one outer fold's best inner-CV score at
% every subset size. Within a fold, gray marks training Delta BIC > 2, black
% marks Delta BIC <= 2, and red marks the eligible subset with the highest
% inner-CV score. Exact score ties keep the smaller subset.

thisDir = fileparts(mfilename('fullpath'));
detailDir = fullfile(thisDir, 'outputs');
figureDir = fullfile(detailDir, 'response_figures');
if ~isfolder(figureDir), mkdir(figureDir); end

threshold = 2;
modelCodes = ["ridge","lasso","elasticnet"];
modelLabels = ["Ridge","Lasso","Elastic net"];
eligibleColor = [0 0 0];
selectColor = [193 18 31] / 255;
excludedColor = [0.78 0.78 0.78];
markerSize = 110;
excludedAlpha = 0.55;
eligibleAlpha = 0.82;
selectAlpha = 0.55;

fig = figure('Color', 'w', 'Units', 'pixels', ...
    'Position', [40 40 1750 520], 'Visible', 'off', ...
    'MenuBar', 'none', 'ToolBar', 'none');
tl = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
rowAxes = gobjects(1, 3);
rowLim = [inf -inf];

for m = 1:numel(modelCodes)
    stem = "mrdflow__markov__" + modelCodes(m);
    detailFile = fullfile(detailDir, char(stem + "__outer_fold_details.csv"));
    D = readtable(detailFile, 'TextType', 'string', 'VariableNamingRule', 'preserve');
    S = chooseWithinOuterFolds(D, threshold);
    writetable(S, fullfile(figureDir, ...
        char(stem + "__deltaBIC2__fold_selections.csv")));
    ks = unique(D.k);

    ax = nexttile(tl);
    rowAxes(m) = ax;
    hold(ax, 'on');
    R = cell(height(S), 1);
    for f = 1:height(S)
        R{f} = sortrows(D(D.repeat == S.repeat(f) & D.outerFold == S.outerFold(f), :), 'k');
        plot(ax, R{f}.k, R{f}.innerSelectionValue, '-', ...
            'Color', [0.45 0.45 0.45 0.45], 'LineWidth', 1.0);
    end
    for f = 1:height(S)
        excluded = R{f}.trainingBIC - min(R{f}.trainingBIC) > threshold + 1e-12;
        scatter(ax, R{f}.k(excluded), R{f}.innerSelectionValue(excluded), ...
            markerSize, excludedColor, 'filled', ...
            'MarkerFaceAlpha', excludedAlpha, 'MarkerEdgeColor', 'none');
    end
    for f = 1:height(S)
        excluded = R{f}.trainingBIC - min(R{f}.trainingBIC) > threshold + 1e-12;
        keep = ~excluded & R{f}.k ~= S.k(f);
        scatter(ax, R{f}.k(keep), R{f}.innerSelectionValue(keep), ...
            markerSize, eligibleColor, 'filled', ...
            'MarkerFaceAlpha', eligibleAlpha, 'MarkerEdgeColor', 'none');
    end
    scatter(ax, S.k, S.innerSelectionValue, markerSize, selectColor, 'filled', ...
        'MarkerFaceAlpha', selectAlpha, 'MarkerEdgeColor', 'none');
    rowLim = [min(rowLim(1), min(D.innerSelectionValue)), ...
        max(rowLim(2), max(D.innerSelectionValue))];

    title(ax, char(modelLabels(m)), 'FontWeight', 'bold', 'FontSize', 14, 'Color', 'k');
    if m == 1
        ylabel(ax, {'Markov', 'Inner-CV score'});
    end
    xlabel(ax, 'Number of mechanistic features');
    ax.XLim = [0.5 max(ks) + 0.5];
    ax.XTick = ks;
    ax.TickLength = [0 0];
    ax.FontName = 'Arial';
    ax.FontSize = 14;
    ax.FontWeight = 'bold';
    ax.XTickLabelRotation = 0;
    ax.Box = 'on';
    ax.YGrid = 'on';
    ax.GridAlpha = 0.14;
end

pad = 0.04 * diff(rowLim);
set(rowAxes, 'YLim', [rowLim(1) - pad, rowLim(2) + pad]);
title(tl, 'Balanced MRD VP', 'FontName', 'Arial', 'FontSize', 16, ...
    'FontWeight', 'bold', 'Interpreter', 'none');

axsAll = findall(fig, 'Type', 'axes');
for i = 1:numel(axsAll)
    if isprop(axsAll(i), 'Toolbar') && ~isempty(axsAll(i).Toolbar)
        axsAll(i).Toolbar.Visible = 'off';
    end
    disableDefaultInteractivity(axsAll(i));
end

stemOut = 'Figure_balanced_MRD_VP_fold_overlay';
exportgraphics(fig, fullfile(figureDir, [stemOut '.png']), 'Resolution', 300);
exportgraphics(fig, fullfile(figureDir, [stemOut '.pdf']), 'ContentType', 'vector');
close(fig);
fprintf('Wrote %s\n', fullfile(figureDir, [stemOut '.png']));
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
