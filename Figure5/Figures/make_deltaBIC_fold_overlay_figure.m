function make_deltaBIC_fold_overlay_figure()
%MAKE_DELTABIC_FOLD_OVERLAY_FIGURE
% One panel per outcome x model with the 10 outer folds drawn on top of
% each other. Each line is one fold's best inner-CV score at every k.
%   gray dot : Delta BIC > 2 in that fold (excluded)
%   black    : Delta BIC <= 2 (eligible)
%   red dot  : highest eligible score in the fold (ties -> smaller k)
% One figure per clinical status. Rows: Flow, Markov.
% Columns: Ridge, Lasso, Elastic net.

thisDir = fileparts(mfilename('fullpath'));
figureRoot = fileparts(thisDir);
addpath(figureRoot);
detailDirs = {
    fullfile(figureRoot, 'outputs')
    fullfile(figureRoot, 'outputs_10fold')
    };
selectionDirs = {
    fullfile(figureRoot, 'outputs', 'delta_bic_selection')
    fullfile(figureRoot, 'outputs_10fold', 'delta_bic_selection')
    };
figureDir = fullfile(figureRoot, 'outputs', 'figures');
if ~isfolder(figureDir), mkdir(figureDir); end
threshold = 2;
views = {
    "flow", "Flow proportion", "Number of flow proportions"
    "markov", "Markov", "Number of mechanistic features"
    };

modelCodes = ["ridge", "lasso", "elasticnet"];
modelLabels = ["Ridge", "Lasso", "Elastic net"];
eligibleColor = [0 0 0];
jobs = {
    "bcr", "BCR::ABL1"
    "mrdflow", "MRD"
    };
selectColor = [193 18 31] / 255;
excludedColor = [0.78 0.78 0.78];
markerSize = 110;
excludedAlpha = 0.55;
eligibleAlpha = 0.82;
selectAlpha = 0.55;

for o = 1:size(jobs, 1)
    outcomeCode = jobs{o, 1};
    outcomeLabel = jobs{o, 2};
    missing = strings(0, 1);
    for v = 1:size(views, 1)
        for m = 1:numel(modelCodes)
            stem = outcomeCode + "__" + views{v, 1} + "__" + modelCodes(m);
            if findNamedFile(detailDirs, stem + "__outer_fold_details.csv") == ""
                missing(end + 1) = stem + "__outer_fold_details.csv"; %#ok<AGROW>
            end
            if findNamedFile(selectionDirs, stem + "__deltaBIC2__fold_selections.csv") == ""
                missing(end + 1) = stem + "__deltaBIC2__fold_selections.csv"; %#ok<AGROW>
            end
        end
    end
    if ~isempty(missing)
        fprintf('Skipping %s overlay; missing saved tables:\n', outcomeLabel);
        fprintf('  %s\n', missing);
        continue
    end
fig = figure('Color', 'w', 'Units', 'pixels', ...
    'Position', [40 40 1750 920], 'Visible', 'off', ...
    'MenuBar', 'none', 'ToolBar', 'none');
tl = tiledlayout(fig, size(views, 1), 3, 'TileSpacing', 'compact', 'Padding', 'compact');

for v = 1:size(views, 1)
    viewCode = views{v, 1};
    viewLabel = views{v, 2};
    xLabel = views{v, 3};
    rowAxes = gobjects(1, 3);
    rowLim = [inf -inf];
    for m = 1:numel(modelCodes)
        stem = outcomeCode + "__" + viewCode + "__" + modelCodes(m);
        D = readtable(findNamedFile(detailDirs, stem + "__outer_fold_details.csv"), ...
            'TextType', 'string', 'VariableNamingRule', 'preserve');
        S = readtable(findNamedFile(selectionDirs, ...
            stem + "__deltaBIC2__fold_selections.csv"), ...
            'TextType', 'string', 'VariableNamingRule', 'preserve');
        ks = unique(D.k);

        ax = nexttile(tl);
        rowAxes(m) = ax;
        hold(ax, 'on');
        % lines first, then dots by layer so red picks stay on top
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

        if v == 1
            title(ax, char(modelLabels(m)), 'FontWeight', 'bold', ...
                'FontSize', 14, 'Color', 'k');
        end
        if m == 1
            ylabel(ax, {char(viewLabel), 'Inner-CV score'});
        end
        xlabel(ax, xLabel);
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
end

title(tl, outcomeLabel, 'FontName', 'Arial', 'FontSize', 16, ...
    'FontWeight', 'bold', 'Interpreter', 'none');

axsAll = findall(fig, 'Type', 'axes');
for i = 1:numel(axsAll)
    if isprop(axsAll(i), 'Toolbar') && ~isempty(axsAll(i).Toolbar)
        axsAll(i).Toolbar.Visible = 'off';
    end
    disableDefaultInteractivity(axsAll(i));
end

if outcomeCode == "bcr"
    stemOut = 'Figure_deltaBIC_BCR_fold_overlay';
else
    stemOut = 'Figure_deltaBIC_MRD_fold_overlay';
end
exportgraphics(fig, fullfile(figureDir, [stemOut '.png']), 'Resolution', 300);
close(fig);
fprintf('Wrote %s\n', fullfile(figureDir, [stemOut '.png']));
end
end

function found = findNamedFile(searchDirs, name)
found = "";
fileName = char(name);
for i = 1:numel(searchDirs)
    candidate = fullfile(searchDirs{i}, fileName);
    if exist(candidate, 'file') == 2
        found = string(candidate);
        return
    end
end
end
