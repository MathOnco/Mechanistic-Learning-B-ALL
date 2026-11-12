function make_response_supp_figures(whichFigure)
% Supplementary figures for the revision response.
% Figure S_BCR_flowMRD: BCR::ABL1 versus flow-MRD among diagnosis specimens.
% Figure S_T4_T5_T6: mean Markov rates from supplemental tables T4, T5, and T6.

    if nargin < 1
        whichFigure = "all";
    end
    outDir = fileparts(mfilename("fullpath"));
    set(0, "DefaultAxesFontName", "Arial");
    set(0, "DefaultTextFontName", "Arial");

    if whichFigure == "all" || whichFigure == "fisher"
        makeFisherFigure(outDir);
    end
    if whichFigure == "all" || whichFigure == "transitions"
        makeTransitionFigure(outDir);
    end
end

function makeFisherFigure(outDir)
    % Same eligibility as the classifiers.
    % BCR::ABL1 cohort: diagnosis, BCR::ABL1 or Negative, finite predictors (n = 101).
    % Flow-MRD cohort: diagnosis, a flow-MRD result, finite predictors (n = 104).
    % Flow-MRD positive = Positive (> 10^-3). Positive (< 10^-3) is negative.
    % Five BCR specimens have no flow-MRD result. Eight flow-MRD specimens are BCR::ABL1-like.
    fisherCounts = [7 42; 9 38];
    [~, p] = fishertest(fisherCounts, "Tail", "both");
    cohortCounts = [7 42; 9 38; 3 5];

    fig = figure("Color", "w", "Units", "pixels", "Position", [80 80 1760 860], "Visible", "off");
    tiledlayout(1, 2, "Padding", "compact", "TileSpacing", "compact");

    nexttile;
    withResult = [49 47];
    positive = [7 9];
    missing = [1 4];
    pct = 100 * positive ./ withResult;
    b = bar(pct, 0.62, "FaceColor", "flat");
    b.CData = [65 166 89; 51 115 191] / 255;
    set(gca, "XTickLabel", {"BCR-negative", "BCR-positive"}, ...
        "FontSize", 18, "FontWeight", "bold", "LineWidth", 1.2, "Box", "off");
    xlabel("BCR::ABL1 cohort (n = 101)", "FontSize", 18, "FontWeight", "bold");
    ylabel("Flow-MRD positive (%)", "FontSize", 18, "FontWeight", "bold");
    ylim([0 40]);
    yticks(0:10:40);
    callouts = {"7/49", "9/47"};
    for i = 1:2
        text(i, pct(i) + 1.6, sprintf("%.1f%%", pct(i)), ...
            "HorizontalAlignment", "center", "FontSize", 20, "FontWeight", "bold");
        text(i, pct(i) + 4.4, callouts{i}, ...
            "HorizontalAlignment", "center", "FontSize", 16, "FontWeight", "bold");
    end
    text(1.5, 36, "No flow-MRD result: 1 negative, 4 positive", ...
        "HorizontalAlignment", "center", "FontSize", 14, "FontWeight", "bold");
    title("50 negative, 51 positive", "FontSize", 22, "FontWeight", "bold");

    nexttile;
    axis([0 2 0 3]);
    axis ij;
    hold on;
    for r = 0:2
        for c = 0:1
            rectangle("Position", [c, r, 1, 1], "FaceColor", [0.96 0.97 0.99], ...
                "EdgeColor", [0.2 0.2 0.2], "LineWidth", 1.2);
            text(c + 0.5, r + 0.5, sprintf("%d", cohortCounts(r + 1, c + 1)), ...
                "HorizontalAlignment", "center", "FontSize", 28, "FontWeight", "bold");
        end
    end
    hold off;
    set(gca, "XTick", [0.5 1.5], "YTick", [0.5 1.5 2.5], ...
        "XTickLabel", {"Flow-MRD positive", "Flow-MRD negative"}, ...
        "YTickLabel", {"BCR-negative", "BCR-positive", "BCR::ABL1-like"}, ...
        "FontSize", 15, "FontWeight", "bold", "TickLength", [0 0], "Box", "off");
    xlabel("Flow-MRD cohort (n = 104)", "FontSize", 18, "FontWeight", "bold");
    title(sprintf("Fisher's exact test, p = %.2f", p), "FontSize", 22, "FontWeight", "bold");

    exportgraphics(fig, fullfile(outDir, "FigureS_BCR_flowMRD.png"), "Resolution", 300);
    close(fig);
    fprintf("Fisher p = %.6f\n", p);
end

function makeTransitionFigure(outDir)
    features = ["M11","M12","M13","M14","M21","M22","M23","M24", ...
                "M31","M32","M33","M34","M41","M42","M43","M44"];
    remission = [0.05 0.36 0.55 0.04 0.02 0.38 0.59 0.02 0.02 0.37 0.60 0.01 0.05 0.35 0.55 0.05];
    remissionSd = [0.01 0.16 0.16 0.02 0.02 0.18 0.18 0.02 0.02 0.18 0.19 0.01 0.01 0.16 0.16 0.02];
    relapse = [0.27 0.35 0.32 0.06 0.23 0.38 0.34 0.05 0.22 0.36 0.37 0.05 0.24 0.35 0.33 0.08];
    relapseSd = [0.34 0.25 0.29 0.10 0.31 0.28 0.31 0.10 0.29 0.25 0.31 0.09 0.28 0.23 0.28 0.09];
    bm = [0.27 0.37 0.31 0.05 0.22 0.41 0.33 0.04 0.22 0.38 0.36 0.04 0.24 0.37 0.33 0.07];
    bmSd = [0.33 0.26 0.29 0.07 0.31 0.28 0.31 0.07 0.29 0.25 0.31 0.07 0.28 0.24 0.28 0.07];
    pb = [0.57 0.25 0.13 0.05 0.53 0.28 0.14 0.05 0.49 0.26 0.17 0.07 0.51 0.26 0.15 0.08];
    pbSd = [0.36 0.24 0.25 0.10 0.33 0.23 0.26 0.10 0.31 0.21 0.27 0.10 0.30 0.21 0.24 0.09];

    fig = figure("Color", "w", "Units", "pixels", "Position", [40 40 1900 820], "Visible", "off");
    tiledlayout(1, 2, "Padding", "compact", "TileSpacing", "compact");

    nexttile;
    drawGrouped(features, remission, remissionSd, relapse, relapseSd, ...
        "Remission (n = 38)", "Relapse (n = 59)", [0.75 0.22 0.55], [51 115 191] / 255);
    title("Summary statistics of remission and relapse", "FontSize", 20, "FontWeight", "bold");

    nexttile;
    drawGrouped(features, bm, bmSd, pb, pbSd, ...
        "Diagnosis BM (n = 109)", "Diagnosis PB (n = 46)", [65 166 89] / 255, [0.80 0.20 0.20]);
    title("Summary statistics of bone marrow and peripheral blood", "FontSize", 20, "FontWeight", "bold");

    exportgraphics(fig, fullfile(outDir, "FigureS_T4_T5_T6.png"), "Resolution", 300);
    close(fig);
end

function drawGrouped(features, a, aSd, b, bSd, nameA, nameB, colorA, colorB)
    x = 1:numel(features);
    hold on;
    eb1 = errorbar(x - 0.16, a, aSd, "o", "LineStyle", "none", ...
        "MarkerFaceColor", colorA, "MarkerEdgeColor", colorA, "Color", colorA, ...
        "LineWidth", 1.4, "MarkerSize", 7, "CapSize", 6);
    eb2 = errorbar(x + 0.16, b, bSd, "o", "LineStyle", "none", ...
        "MarkerFaceColor", colorB, "MarkerEdgeColor", colorB, "Color", colorB, ...
        "LineWidth", 1.4, "MarkerSize", 7, "CapSize", 6);
    hold off;
    xlim([0.4 numel(features) + 0.6]);
    ylim([0 1.35]);
    yticks(0:0.2:1);
    xticks(x);
    xticklabels(features);
    set(gca, "FontSize", 14, "FontWeight", "bold", "LineWidth", 1.2, "Box", "off", ...
        "XTickLabelRotation", 45);
    ylabel("Mean transition rate", "FontSize", 18, "FontWeight", "bold");
    legend([eb1, eb2], {nameA, nameB}, "Location", "northeast", "FontSize", 14, "Box", "off");
end
