function histogramPatchVelocities(histo_v_cell, name, locsLabel, M0_ff_img)

ToolBox = getGlobalToolBox;

% Check sizes
[rows, cols] = size(locsLabel);
[numX, numY] = size(M0_ff_img);
assert(isequal(size(histo_v_cell), [rows, cols]), 'Size of histo_v_cell must match locsLabel');

% Create figure and display image
fi = figure("Visible", "on", 'Color', 'w', ...
    'Units', 'pixels', "Position", [200 200 600 600]);

% Create main axes for image with controlled position
ax_main = axes('Position', [0.1 0.1 0.8 0.8]); % Leave 10% margin on all sides
imshow(M0_ff_img, [], 'Parent', ax_main);
axis(ax_main, 'image');
axis(ax_main, 'off');
title(ax_main, ['Velocity Histograms Overlay - ' name]);

% Parameters for histogram size (in normalized units relative to image)
histWidth = 0.04;  % width (5% of image width)
histHeight = 0.03;  % height (4% of image height)

% Get axis limits (important for coordinate conversion)
xlim = ax_main.XLim;
ylim = ax_main.YLim;

% Calculate conversion factors
x_scale = diff(xlim)/numX;
y_scale = diff(ylim)/numY;

% Store the axis position in normalized figure units
ax_pos = ax_main.Position; % [left, bottom, width, height]

for circleIdx = 1:rows

    for i = 1:cols

        if isempty(locsLabel{circleIdx, i}) || isempty(histo_v_cell{circleIdx, i})
            continue;
        end

        % Get histogram data
        histData = histo_v_cell{circleIdx, i};
        if numel(histData) ~= 2
            continue;
        end

        counts = histData{1};
        edges = histData{2};
        if isempty(counts) || isempty(edges)
            continue;
        end

        % Compute histogram center location
        pos = locsLabel{circleIdx, i}; % pos = [x, y]
        if isempty(pos) || numel(pos) ~= 2
            continue;
        end

        % Convert image coordinates to axes coordinates
        x_img = xlim(1) + (pos(1)-1)*x_scale;  % -1 because MATLAB is 1-indexed
        y_img = ylim(2) - (pos(2)-1)*y_scale;  % Flip y-axis

        % Create histogram axes (centered on point)
        ax_hist = axes('Position', ...
            [ax_main.Position(1) + (x_img-xlim(1))/diff(xlim)*ax_main.Position(3) - histWidth/2, ...
            ax_main.Position(2) + (y_img-ylim(1))/diff(ylim)*ax_main.Position(4) - histHeight/2, ...
            histWidth, ...
            histHeight], ...
            'Units', 'normalized');
        % Plot histogram
        bar(ax_hist, edges(1:end - 1), counts, 'histc');
        ax_hist.XTick = [];
        ax_hist.YTick = [];
        ax_hist.Box = 'off';
        ax_hist.Color = 'none';
        ax_hist.XColor = 'none';
        ax_hist.YColor = 'none';
    end

end

outputDir = fullfile(ToolBox.path_png, 'local');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

% Save figure
exportgraphics(gcf, fullfile(outputDir, ...
    sprintf("%s_velocities_histogram_overlay_%s.png", ToolBox.folder_name, name)));

close(fi);

end
