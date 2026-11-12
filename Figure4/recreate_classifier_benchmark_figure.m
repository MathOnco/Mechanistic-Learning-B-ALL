function files = recreate_classifier_benchmark_figure(csvFile, outputDir)
%RECREATE_CLASSIFIER_BENCHMARK_FIGURE Recreate the nested-CV CI figure.
%
% files = recreate_classifier_benchmark_figure()
% files = recreate_classifier_benchmark_figure(csvFile, outputDir)
%
% The input CSV must contain these columns:
%   analysis, model, metric, estimate, ciLower, ciUpper
%
% The expected analysis labels have the form OUTCOME__VIEW, for example
% BCR__markov or MRDflow__flow. The script plots AUC and balanced accuracy
% for the four predictor views and four classifiers used in the benchmark.

thisDir = fileparts(mfilename('fullpath'));
if nargin < 1 || isempty(csvFile)
    csvFile = fullfile(thisDir, 'benchmark_all_bootstrap_95CI.csv');
end
if nargin < 2 || isempty(outputDir)
    outputDir = fullfile(thisDir, 'outputs');
end
csvFile = char(string(csvFile));
outputDir = char(string(outputDir));

if ~isfile(csvFile)
    error('benchmarkFigure:MissingCSV', ...
        'Input CSV was not found: %s', csvFile);
end
if ~isfolder(outputDir)
    mkdir(outputDir);
end

T = readtable(csvFile, 'TextType', 'string', ...
    'VariableNamingRule', 'preserve');
required = ["analysis","model","metric","estimate","ciLower","ciUpper"];
missing = required(~ismember(required, string(T.Properties.VariableNames)));
if ~isempty(missing)
    error('benchmarkFigure:MissingColumns', ...
        'CSV is missing required column(s): %s', strjoin(missing, ', '));
end

T.analysis = string(T.analysis);
T.model = lower(string(T.model));
T.metric = string(T.metric);

% The attached layout places BCR on the first row and Flow MRD second.
outcomes = ["BCR", "MRDflow"];
metrics = ["auc", "balancedAccuracy"];
views = ["flow", "markov", "clinical_genomic", "markov_plus_clinical"];
models = ["ridge", "lasso", "elasticnet", "randomforest"];

fig = figure('Color', 'w', 'Units', 'pixels', ...
    'Position', [60 60 1500 920], 'Visible', 'off');
tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', ...
    'Padding', 'compact');

for o = 1:numel(outcomes)
    for q = 1:numel(metrics)
        ax = nexttile(tl, (o - 1) * 2 + q);
        [S, labels, colors] = orderedPanelRows( ...
            T, outcomes(o), metrics(q), views, models);
        drawPanel(ax, S, labels, colors);
    end
end

pngFile = fullfile(outputDir, 'Figure_R1_classifier_benchmark.png');
pdfFile = fullfile(outputDir, 'Figure_R1_classifier_benchmark.pdf');
exportgraphics(fig, pngFile, 'Resolution', 300);
exportgraphics(fig, pdfFile, 'ContentType', 'vector');
close(fig);

files = struct('png', string(pngFile), 'pdf', string(pdfFile), ...
    'sourceCSV', string(csvFile));
fprintf('Created:\n  %s\n  %s\n', pngFile, pdfFile);
end

function [S, labels, colors] = orderedPanelRows(T, outcome, metric, views, models)
% Enforce identical row ordering across all four panels.
nExpected = numel(views) * numel(models);
S = T([], :);
labels = strings(nExpected, 1);
colors = zeros(nExpected, 3);
row = 0;

for v = 1:numel(views)
    for m = 1:numel(models)
        row = row + 1;
        analysisName = outcome + "__" + views(v);
        hit = T.analysis == analysisName & T.model == models(m) & ...
            T.metric == metric;
        if sum(hit) ~= 1
            error('benchmarkFigure:UnexpectedRowCount', ...
                ['Expected one row for analysis=%s, model=%s, metric=%s; ' ...
                 'found %d.'], analysisName, models(m), metric, sum(hit));
        end
        S = [S; T(hit, :)]; %#ok<AGROW>
        labels(row) = prettyView(views(v)) + " | " + prettyModel(models(m));
        colors(row, :) = viewColor(views(v));
    end
end
end

function drawPanel(ax, S, labels, colors)
n = height(S);
y = (n:-1:1)';
hold(ax, 'on');

for i = 1:n
    lowerError = S.estimate(i) - S.ciLower(i);
    upperError = S.ciUpper(i) - S.estimate(i);
    errorbar(ax, S.estimate(i), y(i), lowerError, upperError, ...
        'horizontal', 'o', 'Color', colors(i, :), ...
        'MarkerFaceColor', colors(i, :), 'MarkerEdgeColor', colors(i, :), ...
        'LineWidth', 1.4, 'MarkerSize', 5, 'CapSize', 0);
end

xline(ax, 0.5, ':', 'Chance', 'Color', [0.45 0.45 0.45], ...
    'LabelOrientation', 'aligned', 'LabelVerticalAlignment', 'bottom', ...
    'FontSize', 13);
xlim(ax, [0 1]);
ylim(ax, [0.5 n + 0.5]);
yticks(ax, 1:n);
yticklabels(ax, flipud(labels));
set(ax, 'YDir', 'normal', 'FontSize', 13, 'TickLabelInterpreter', 'none', ...
    'Box', 'off', 'Layer', 'top');
grid(ax, 'on');
ax.GridAlpha = 0.16;
ax.MinorGridAlpha = 0.08;
hold(ax, 'off');
end

function s = prettyView(s)
switch s
    case "markov"
        s = "Markov";
    case "flow"
        s = "Flow proportion";
    case "clinical_genomic"
        s = "Clinical/genomic";
    case "markov_plus_clinical"
        s = "Markov + clinical/genomic";
end
end

function s = prettyModel(s)
switch s
    case "ridge"
        s = "Ridge";
    case "lasso"
        s = "Lasso";
    case "elasticnet"
        s = "Elastic net";
    case "randomforest"
        s = "Random forest";
end
end

function c = viewColor(view)
switch view
    case "markov"
        c = [0.20 0.45 0.75];
    case "flow"
        c = [0.25 0.65 0.35];
    case "clinical_genomic"
        c = [0.85 0.45 0.15];
    otherwise
        c = [0.55 0.30 0.70];
end
end
