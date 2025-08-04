function [f] = velocityIm(v_mean, mask, cmap, name, options)

arguments
    v_mean
    mask
    cmap
    name
    options.colorbarOn
    options.LabelName = 'mm/s'
end

ToolBox = getGlobalToolBox;

f = figure("Visible", "off");
colormap(cmap)
imagesc(v_mean .* mask)

if options.colorbarOn
    c = colorbar;
    c.Label.String = options.LabelName;
end

axis image; axis off;

% Export the figure as a PNG with a white background
exportgraphics(gcf, fullfile(ToolBox.path_png, sprintf("%s_map_%s.png", ToolBox.folder_name, name)), 'Resolution', 300);

end
