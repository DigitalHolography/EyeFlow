function fig = figAspect(opt)

arguments
    opt.FigName = ''
    opt.Visible = 'off'
    opt.Fontsize = 12
    opt.LineWidth = 2
end

fig = figure('Name', opt.FigName, 'Visible', opt.Visible, 'Color', 'w');

fontsize(gca, opt.Fontsize, "points");

pbaspect([1.618 1 1]);

box on
set(gca, 'LineWidth', opt.LineWidth);

axis tight;

end
