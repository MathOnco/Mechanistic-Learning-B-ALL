function exportResponseFigure(fig, P, stem)
%EXPORTRESPONSEFIGURE Export high-resolution PNG and vector PDF.

set(fig, 'Color', 'w');
if P.isSmoke
    annotation(fig, 'textbox', [0.79 0.965 0.19 0.025], ...
        'String', 'SMOKE TEST — NOT REPORTABLE', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
        'Color', [0.75 0.1 0.1], 'FontWeight', 'bold', 'FontSize', 8, ...
        'EdgeColor', 'none', 'BackgroundColor', [1 0.94 0.94]);
end
exportgraphics(fig, fullfile(P.figureDir, [stem '.png']), 'Resolution', 300);
exportgraphics(fig, fullfile(P.figureDir, [stem '.pdf']), 'ContentType', 'vector');
end
