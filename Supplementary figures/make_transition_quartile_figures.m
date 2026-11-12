function make_transition_quartile_figures()
% Quartile plots of all 16 Markov transitions by BCR::ABL1 and by flow-MRD.
% Writes original-specimen and virtual-patient versions (PNG only).

    outDir = fileparts(mfilename("fullpath"));
    dataDir = fullfile(outDir, "..", "Figure4", "Data");
    vpFile = fullfile(outDir, "..", "Figure6", "Data", "PB_BM_merged_with_MRDflow_BCR.csv");
    set(0, "DefaultAxesFontName", "Arial");
    set(0, "DefaultTextFontName", "Arial");

    features = ["M11","M12","M13","M14","M21","M22","M23","M24", ...
                "M31","M32","M33","M34","M41","M42","M43","M44"];

    colorBcrNeg = [235 65 38] / 255;
    colorBcrPos = [16 49 244] / 255;
    colorMrdNeg = [1 1 1];
    colorMrdPos = [47 224 131] / 255;

    [Xbcr0, yBcr0] = loadOriginalOutcomeMatrix(dataDir, features, "bcr");
    [Xmrd0, yMrd0] = loadOriginalOutcomeMatrix(dataDir, features, "mrdflow");
    [XbcrV, yBcrV] = loadVpOutcomeMatrix(vpFile, features, "bcr");
    [XmrdV, yMrdV] = loadVpOutcomeMatrix(vpFile, features, "mrdflow");

    plotQuartiles(Xbcr0, yBcr0, features, ...
        "BCR-negative", "BCR::ABL1-positive", colorBcrNeg, colorBcrPos, ...
        "Markov transition rates by BCR::ABL1 status", ...
        fullfile(outDir, "FigureS_transitions_by_BCR_original.png"));

    plotQuartiles(Xmrd0, yMrd0, features, ...
        "Flow-MRD negative", "Flow-MRD positive", colorMrdNeg, colorMrdPos, ...
        "Markov transition rates by flow-MRD status", ...
        fullfile(outDir, "FigureS_transitions_by_MRD_original.png"), true);

    plotQuartiles(XbcrV, yBcrV, features, ...
        "BCR-negative", "BCR::ABL1-positive", colorBcrNeg, colorBcrPos, ...
        "Virtual-patient Markov transition rates by BCR::ABL1 status", ...
        fullfile(outDir, "FigureS_transitions_by_BCR.png"));

    plotQuartiles(XmrdV, yMrdV, features, ...
        "Flow-MRD negative", "Flow-MRD positive", colorMrdNeg, colorMrdPos, ...
        "Virtual-patient Markov transition rates by flow-MRD status", ...
        fullfile(outDir, "FigureS_transitions_by_MRD.png"), true);
end

function [X, y] = loadOriginalOutcomeMatrix(dataDir, features, outcome)
    Dm = readTransposed(fullfile(dataDir, "markov.csv"));
    Df = readTransposed(fullfile(dataDir, "flow.csv"));
    if Dm.nSamples ~= Df.nSamples
        error("markov.csv and flow.csv have different sample counts.");
    end

    Xall = nan(Dm.nSamples, numel(features));
    for j = 1:numel(features)
        Xall(:, j) = getNumericRow(Dm, features(j))';
    end
    flowNames = ["x1","x2","x3","x4"];
    Xflow = nan(Df.nSamples, numel(flowNames));
    for j = 1:numel(flowNames)
        Xflow(:, j) = getNumericRow(Df, flowNames(j))';
    end

    progression = getStringRow(Dm, "DiseaseProgression")';
    bcr = getStringRow(Dm, "BCR::ABL1/like")';
    mrdFlow = getStringRow(Dm, "MRD(flow)")';
    isDiagnosis = ~strcmpi(progression, "Remission") & ~strcmpi(progression, "Relapse");
    finitePredictors = all(isfinite(Xall), 2) & all(isfinite(Xflow), 2);

    switch lower(string(outcome))
        case "bcr"
            eligible = bcr == "BCR::ABL1" | bcr == "Negative";
            yAll = double(bcr == "BCR::ABL1");
        case "mrdflow"
            eligible = mrdFlow == "Positive (> 10^-3)" | mrdFlow == "Negative" | ...
                mrdFlow == "Positive (< 10^-3)";
            yAll = double(mrdFlow == "Positive (> 10^-3)");
        otherwise
            error("Unknown outcome: %s", outcome);
    end

    keep = isDiagnosis & eligible & finitePredictors;
    X = Xall(keep, :);
    y = yAll(keep);
    fprintf("%s original: n = %d (positive = %d, negative = %d)\n", ...
        outcome, numel(y), sum(y == 1), sum(y == 0));
end

function [X, y] = loadVpOutcomeMatrix(vpFile, features, outcome)
    T = readtable(vpFile, "TextType", "string", "VariableNamingRule", "preserve");
    required = [features, "x1", "x2", "x3", "x4", "MRD(flow)", "BCR::ABL1/like"];
    missing = setdiff(required, string(T.Properties.VariableNames));
    if ~isempty(missing)
        error("Missing columns in VP table: %s", strjoin(missing, ", "));
    end

    Xall = nan(height(T), numel(features));
    for j = 1:numel(features)
        Xall(:, j) = double(T{:, features(j)});
    end
    Xflow = [double(T.x1), double(T.x2), double(T.x3), double(T.x4)];
    finitePredictors = all(isfinite(Xall), 2) & all(isfinite(Xflow), 2);

    bcr = strtrim(string(T{:, "BCR::ABL1/like"}));
    mrd = strtrim(string(T{:, "MRD(flow)"}));

    switch lower(string(outcome))
        case "bcr"
            eligible = ismember(bcr, ["BCR::ABL1", "Negative"]);
            yAll = double(bcr == "BCR::ABL1");
        case "mrdflow"
            eligible = ismember(mrd, ["Positive (> 10^-3)", "Negative", "Positive (< 10^-3)"]);
            yAll = double(mrd == "Positive (> 10^-3)");
        otherwise
            error("Unknown outcome: %s", outcome);
    end

    keep = eligible & finitePredictors;
    X = Xall(keep, :);
    y = yAll(keep);
    fprintf("%s VP: n = %d (positive = %d, negative = %d)\n", ...
        outcome, numel(y), sum(y == 1), sum(y == 0));
end

function plotQuartiles(X, y, features, nameNeg, namePos, colorNeg, colorPos, figTitle, outPath, outlineNeg)
    if nargin < 10
        outlineNeg = false;
    end
    nFeat = numel(features);
    values = X(:);
    featureIdx = repelem((1:nFeat)', size(X, 1));
    group = repmat(categorical(y, [0 1], {char(nameNeg), char(namePos)}), nFeat, 1);

    fig = figure("Color", "w", "Units", "pixels", "Position", [40 40 2000 900], "Visible", "off");
    ax = axes(fig);
    b = boxchart(ax, featureIdx, values, "GroupByColor", group, ...
        "BoxWidth", 0.72, "MarkerStyle", ".", "MarkerSize", 8, "LineWidth", 1.3);

    colororder(ax, [colorNeg; colorPos]);
    if outlineNeg
        b(1).BoxFaceColor = colorNeg;
        b(1).BoxEdgeColor = [0.25 0.25 0.25];
        b(1).WhiskerLineColor = [0.25 0.25 0.25];
        b(1).MarkerColor = [0.45 0.45 0.45];
        b(2).BoxFaceColor = colorPos;
        b(2).BoxEdgeColor = colorPos * 0.75;
        b(2).WhiskerLineColor = colorPos * 0.75;
        b(2).MarkerColor = colorPos;
    end
    xlim(ax, [0.4 nFeat + 0.6]);
    ylim(ax, [-0.02 1.05]);
    yticks(ax, 0:0.2:1);
    xticks(ax, 1:nFeat);
    xticklabels(ax, features);
    set(ax, "FontSize", 16, "FontWeight", "bold", "LineWidth", 1.2, "Box", "off", ...
        "XTickLabelRotation", 45);
    ylabel(ax, "Transition rate", "FontSize", 18, "FontWeight", "bold");
    title(ax, figTitle, "FontSize", 22, "FontWeight", "bold");
    legend(ax, {char(nameNeg), char(namePos)}, "Location", "northeast", ...
        "FontSize", 16, "Box", "off");

    exportgraphics(fig, outPath, "Resolution", 300);
    close(fig);
    fprintf("Wrote %s\n", outPath);
end

function D = readTransposed(filePath)
    raw = readcell(filePath);
    labels = string(raw(:, 1));
    values = raw(:, 2:end);
    D.nSamples = size(values, 2);
    D.labels = labels;
    D.values = values;
end

function v = getNumericRow(D, name)
    idx = find(D.labels == string(name), 1);
    if isempty(idx)
        error("Missing numeric row: %s", name);
    end
    raw = D.values(idx, :);
    v = nan(1, numel(raw));
    for k = 1:numel(raw)
        x = raw{k};
        if isnumeric(x) && isscalar(x)
            v(k) = double(x);
        elseif isstring(x) || ischar(x)
            s = strtrim(string(x));
            if s ~= "" && ~strcmpi(s, "NaN")
                v(k) = str2double(s);
            end
        end
    end
end

function v = getStringRow(D, name)
    idx = find(D.labels == string(name), 1);
    if isempty(idx)
        error("Missing string row: %s", name);
    end
    raw = D.values(idx, :);
    v = strings(1, numel(raw));
    for k = 1:numel(raw)
        x = raw{k};
        if isempty(x)
            v(k) = "NaN";
        elseif isnumeric(x)
            if isscalar(x) && isnan(x)
                v(k) = "NaN";
            else
                v(k) = strtrim(string(x));
            end
        else
            s = strtrim(string(x));
            if ismissing(s) || s == "" || strcmpi(s, "NaN")
                v(k) = "NaN";
            else
                v(k) = s;
            end
        end
    end
end
