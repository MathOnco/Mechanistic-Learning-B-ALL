function files = make_vp_outer_oof_bootstrap_figure(varargin)
%MAKE_VP_OUTER_OOF_BOOTSTRAP_FIGURE Original vs VP mean-OOF bootstrap CIs.
%
% Point estimates are metrics of the mean outer out-of-fold probabilities
% across nested-CV repeats. Original intervals resample patients. VP
% intervals are a patient-cluster bootstrap: each resample draws the
% original patients with replacement and keeps all VP rows of each drawn
% patient (2,000 resamples).
%
%   files = make_vp_outer_oof_bootstrap_figure
%   files = make_vp_outer_oof_bootstrap_figure('ComputeVPIfMissing', true)

figureRoot = fileparts(fileparts(mfilename('fullpath')));
ip = inputParser;
addParameter(ip, 'OriginalDir', fullfile(figureRoot, 'Data'));
addParameter(ip, 'VPDir', fullfile(figureRoot, 'outputs', 'specimen_folds'));
addParameter(ip, 'OutputDir', fullfile(figureRoot, 'outputs', 'response_figures'));
addParameter(ip, 'ComputeVPIfMissing', false, @(x) islogical(x) && isscalar(x));
addParameter(ip, 'Visible', 'off', @(x) ismember(string(x), ["on","off"]));
parse(ip, varargin{:});
opt = ip.Results;
if ~isfolder(opt.OutputDir), mkdir(opt.OutputDir); end

models = ["ridge","lasso","elasticnet"];
modelLabels = ["Ridge","Lasso","Elastic net"];
outcomes = ["BCR","MRDflow"];
outcomeLabels = ["BCR::ABL1","MRD"];
metrics = ["balancedAccuracy","auc"];
metricLabels = ["Balanced accuracy (BA)","AUC"];
repOrder = ["Original","VP"];
repColors = [
    51 115 191
    51 115 191
] / 255;

S = table();
for o = 1:numel(outcomes)
    origFile = fullfile(opt.OriginalDir, ...
        char(lower(outcomes(o)) + "__markov_bootstrap_95CI.csv"));
    vpFile = fullfile(opt.VPDir, ...
        char("vp_" + lower(outcomes(o)) + "__markov_cluster_bootstrap_95CI.csv"));
    if ~isfile(vpFile)
        error(['Missing %s\nRun compute_vp_cluster_bootstrap first.'], vpFile);
    end
    assert(isfile(origFile), 'Missing original bootstrap file: %s', origFile);
    Orig = readtable(origFile, 'TextType', 'string', 'VariableNamingRule', 'preserve');
    VP = readtable(vpFile, 'TextType', 'string', 'VariableNamingRule', 'preserve');
    for m = 1:numel(models)
        for q = 1:numel(metrics)
            for r = 1:numel(repOrder)
                if repOrder(r) == "Original"
                    T = Orig;
                    analysis = outcomes(o) + "__markov";
                else
                    T = VP;
                    analysis = "VP__" + outcomes(o) + "__markov";
                end
                hit = T.analysis == analysis & lower(T.model) == models(m) & ...
                    T.metric == metrics(q);
                assert(sum(hit) == 1, 'Expected one row for %s %s %s', ...
                    analysis, models(m), metrics(q));
                R = T(hit, :);
                S = [S; table(outcomes(o), repOrder(r), models(m), metrics(q), ...
                    R.estimate, R.ciLower, R.ciUpper, R.n, ...
                    'VariableNames', {'outcome','representation','model', ...
                    'metric','estimate','ciLower','ciUpper','n'})]; %#ok<AGROW>
            end
        end
    end
end

picked = S;
writetable(picked, fullfile(opt.OutputDir, ...
    'Figure_VP_outer_oof_bootstrap_summary.csv'));

for labeled = [true false]
    for e = 1:numel(outcomes)
        fig = figure('Color', 'w', 'Units', 'pixels', ...
            'Position', [40 40 2100 860], 'Visible', char(opt.Visible), ...
            'MenuBar', 'none', 'ToolBar', 'none');
        tl = tiledlayout(fig, 1, numel(models), 'TileSpacing', 'compact', 'Padding', 'compact');
        for j = 1:numel(models)
            ax = nexttile(tl);
            hold(ax, 'on');
            yTicks = zeros(1, numel(repOrder));
            for v = 1:numel(repOrder)
                yBase = 1 + (v - 1) * 2.70;
                yTicks(v) = yBase + 0.50;
                color = repColors(v, :);
                for t = 1:numel(metrics)
                    hit = picked.outcome == outcomes(e) & ...
                        picked.model == models(j) & ...
                        picked.representation == repOrder(v) & ...
                        picked.metric == metrics(t);
                    R = picked(hit, :);
                    y = yBase + (t - 1) * 1.00;
                    drawInterval(ax, R.estimate, R.ciLower, R.ciUpper, y, color);
                    if labeled
                        text(ax, R.estimate, y + 0.22, ...
                            sprintf('%.2f  [%.2f, %.2f]', ...
                            R.estimate, R.ciLower, R.ciUpper), ...
                            'HorizontalAlignment', 'center', 'FontName', 'Arial', ...
                            'FontSize', 16, 'FontWeight', 'bold', 'Color', 'k');
                        text(ax, R.estimate, y + 0.52, metricLabels(t), ...
                            'HorizontalAlignment', 'center', 'FontName', 'Arial', ...
                            'FontSize', 14, 'FontWeight', 'bold', 'Color', 'k');
                    end
                end
                if labeled
                    nVal = picked.n(picked.outcome == outcomes(e) & ...
                        picked.model == models(j) & picked.representation == repOrder(v) & ...
                        picked.metric == metrics(1));
                    text(ax, 0.355, yBase - 0.42, sprintf('n=%d', nVal), ...
                        'HorizontalAlignment', 'left', 'FontName', 'Arial', ...
                        'FontSize', 15, 'Color', 'k');
                end
            end
            ax.XLim = [0.35 1.00];
            ax.YLim = [0.30 5.85];
            ax.YDir = 'reverse';
            ax.YTick = yTicks;
            ax.YTickLabel = {"Original Markov","VP Markov"};
            ax.XTick = 0.4:0.1:1.0;
            ax.XGrid = 'on';
            ax.YGrid = 'off';
            ax.GridAlpha = 0.16;
            ax.TickLength = [0 0];
            ax.FontName = 'Arial';
            ax.Box = 'on';
            xlabel(ax, {'Nested-CV mean outer-OOF performance', ...
                'Patient bootstrap 95% CI, 2,000 resamples'});
            ax.FontSize = 22;
            ax.FontWeight = 'bold';
            ax.XAxis.FontSize = 22;
            ax.XAxis.FontWeight = 'bold';
            ax.YAxis.FontSize = 18;
            ax.YAxis.FontWeight = 'bold';
            ax.XLabel.FontSize = 13;
            ax.XLabel.FontWeight = 'normal';
            if labeled
                title(ax, modelLabels(j), 'FontWeight', 'bold', 'FontSize', 16, 'Color', 'k');
            end
        end
        if labeled
            title(tl, outcomeLabels(e) + "  |  outer-OOF AUC and balanced accuracy", ...
                'FontName', 'Arial', 'FontSize', 17, 'FontWeight', 'bold', ...
                'Interpreter', 'none');
        end
        axs = findall(fig, 'Type', 'axes');
        for i = 1:numel(axs)
            if isprop(axs(i), 'Toolbar') && ~isempty(axs(i).Toolbar)
                axs(i).Toolbar.Visible = 'off';
            end
            disableDefaultInteractivity(axs(i));
        end
        stem = char("Figure_VP_" + outcomes(e) + "_outer_oof_bootstrap");
        if ~labeled
            stem = [stem '_notext'];
        end
        exportgraphics(fig, fullfile(opt.OutputDir, [stem '.png']), 'Resolution', 300);
        if labeled
            exportgraphics(fig, fullfile(opt.OutputDir, [stem '.pdf']), 'ContentType', 'vector');
        end
        close(fig);
        fprintf('Wrote %s\n', fullfile(opt.OutputDir, [stem '.png']));
    end
end

files.summary = string(fullfile(opt.OutputDir, ...
    'Figure_VP_outer_oof_bootstrap_summary.csv'));
files.bcr = string(fullfile(opt.OutputDir, 'Figure_VP_BCR_outer_oof_bootstrap.png'));
files.mrd = string(fullfile(opt.OutputDir, 'Figure_VP_MRDflow_outer_oof_bootstrap.png'));
end

function drawInterval(ax, point, lower, upper, y, color)
line(ax, [lower upper], [y y], 'Color', color, 'LineWidth', 2.2);
line(ax, [lower lower], [y - 0.08 y + 0.08], 'Color', color, 'LineWidth', 1.6);
line(ax, [upper upper], [y - 0.08 y + 0.08], 'Color', color, 'LineWidth', 1.6);
plot(ax, point, y, 'd', 'MarkerSize', 9, ...
    'MarkerFaceColor', color, 'MarkerEdgeColor', 'w', 'LineWidth', 1.0);
end
