function [f] = graphSignal(filename, folder, x, y, style, color, opt)
% Graph Function for Signal display

arguments
    filename {mustBeText}
    folder {mustBeText}
end

arguments (Repeating)
    x {mustBeNumeric}
    y {mustBeNumeric}
    style
    color
end

arguments
    opt.xlabel {mustBeText} = 'default'
    opt.ylabel {mustBeText} = 'default'
    opt.Title {mustBeText} = 'default'
    opt.LineWidth = 2
    opt.Fontsize = 14
    opt.Legends = {}
    opt.TxtName = {}
    opt.TxtFigX = []
    opt.TxtFigY = []
    opt.TxtFigString = []
    opt.yLines = []
    opt.yLineLabels = {}
    opt.xLines = []
    opt.xLineLabels = {}
end

f = figAspect('Fontsize', opt.Fontsize, 'LineWidth', opt.LineWidth);
hold on

for n = 1:length(y)
    plot(x{n}, y{n}, style{n}, 'Color', color{n}, 'LineWidth', opt.LineWidth)

end

if ~isempty(opt.TxtFigX)

    for n = 1:length(opt.TxtFigX)
        text(opt.TxtFigX(n), opt.TxtFigY(n), opt.TxtFigString{n})
    end

end

if ~isempty(opt.yLines)

    for n = 1:length(opt.yLines)
        yline(opt.yLines(n), ':', opt.yLineLabels{n}, LineWidth = opt.LineWidth);
    end

end

if ~isempty(opt.xLines)

    for n = 1:length(opt.xLines)
        xline(opt.xLines(n), ':', obj.xLineLabels{n}, LineWidth = opt.LineWidth);
    end

end

title(opt.Title)

xlabel(opt.xlabel, 'FontSize', opt.Fontsize);
ylabel(opt.ylabel, 'FontSize', opt.Fontsize);

if ~isempty(opt.Legends)
    legend(opt.Legends)
end

if ~isempty(opt.TxtName)

    for n = 1:length(opt.TxtName)
        plot2txt(x{n}, y{n}, opt.TxtName{n})
    end

end

axis padded
axP = axis;
axis tight
axT = axis;
axis([axT(1), axT(2), axP(3), axP(4)])

ToolBox = getGlobalToolBox;
exportgraphics(gca, fullfile(ToolBox.path_png, folder, sprintf("%s_%s_graph.png", ToolBox.main_foldername, filename)))
exportgraphics(gca, fullfile(ToolBox.path_eps, folder, sprintf("%s_%s_graph.eps", ToolBox.main_foldername, filename)))

end
