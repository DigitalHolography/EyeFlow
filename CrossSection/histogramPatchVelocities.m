function histogramPatchVelocities(histo_v_cell, name, locsLabel, M0_ff_img)

ToolBox = getGlobalToolBox;

params = ToolBox.getParams;
exportVideos = params.exportVideos;

% Check sizes
[rows, cols] = size(locsLabel);
[numX, numY] = size(M0_ff_img);
assert(isequal(size(histo_v_cell), [rows, cols]), 'Size of histo_v_cell must match locsLabel');

i = 0;
numFrames = 0;

while numFrames <= 0
    i = i + 1;
    numFrames = size(histo_v_cell{i}, 2);

    if i > size(histo_v_cell, 1) * size(histo_v_cell, 2)
        warning("Velocity profiles cells are all empty.")
        break
    end

end

% Create figure and display image
fi = figure("Visible", "on", 'Color', 'w', ...
    'Units', 'pixels', "Position", [200 200 600 600]);

% Create main axes for image with controlled position
ax_main = axes('Position', [0.1 0.1 0.8 0.8]); % Leave 10 % margin on all sides
imshow(M0_ff_img, [], 'Parent', ax_main);
axis(ax_main, 'image');
axis(ax_main, 'off');
title(ax_main, ['Velocity Histograms Overlay - ' name]);

% Parameters for histogram size (in normalized units relative to image)
histWidth = 0.04; % width (5 % of image width)
histHeight = 0.03; % height (4 % of image height)

% Get axis limits (important for coordinate conversion)
xlim = ax_main.XLim;
ylim = ax_main.YLim;

% Calculate conversion factors
x_scale = diff(xlim) / numX;
y_scale = diff(ylim) / numY;

% Store the axis position in normalized figure units
ax_pos = ax_main.Position; % [left, bottom, width, height]

plot_list = cell(rows, cols);
for circleIdx = 1:rows

    for i = 1:cols

        if isempty(locsLabel{circleIdx, i}) || isempty(histo_v_cell{circleIdx, i})
            continue;
        end

        % Get histogram data
        histData = zeros(1,5);
        histo_t = histo_v_cell{circleIdx, i};
        for ff = 1:numFrames
            histo = histo_t{ff};
            histData = histData + histo;
        end
        histData = histData / numFrames;

        counts = histData;
        edges = linspace(0,60,6); %% HARD CODED
        if isempty(counts) || isempty(edges)
            continue;
        end

        % Compute histogram center location
        pos = locsLabel{circleIdx, i}; % pos = [x, y]

        if isempty(pos) || numel(pos) ~= 2
            continue;
        end

        % Convert image coordinates to axes coordinates
        x_img = xlim(1) + (pos(1) - 1) * x_scale; % -1 because MATLAB is 1-indexed
        y_img = ylim(2) - (pos(2) - 1) * y_scale; % Flip y-axis

        % Create histogram axes (centered on point)
        ax = axes('Position', ...
            [ax_main.Position(1) + (x_img - xlim(1)) / diff(xlim) * ax_main.Position(3) - histWidth / 2, ...
             ax_main.Position(2) + (y_img - ylim(1)) / diff(ylim) * ax_main.Position(4) - histHeight / 2, ...
             histWidth, ...
             histHeight], ...
            'Units', 'normalized');
        % Plot histogram
        bh = bar(ax, edges(1:end-1), counts,'BarWidth', 1.0);
        plot_list{circleIdx, i} = bh;
        ax.XTick = [];
        ax.YTick = [];
        ax.Box = 'off';
        ax.Color = 'none';
        ax.XColor = 'none';
        ax.YColor = 'none';
    end

end
hold off

outputDir = fullfile(ToolBox.path_png);

if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

% Save figure
saveas(gcf, fullfile(outputDir, ...
    sprintf("%s_histogram_velocities_overlay_%s.png", ToolBox.folder_name, name)));

%% (GIF)
if exportVideos
hold on;

histPatchVelocitiesVideo = zeros(600,600, 3, numFrames,'single');
for frameIdx = 1:numFrames
    %fprintf(" %d ",frameIdx);
    for circleIdx = 1:rows
        for i = 1:cols
            if isempty(locsLabel{circleIdx, i}) || isempty(histo_v_cell{circleIdx, i})
                continue;
            end

            % Get histogram data
            histo_t = histo_v_cell{circleIdx, i};
            histData = histo_t{frameIdx};

            % replot
            % Ensure histData matches number of bars
            bh = plot_list{circleIdx, i};
            if numel(bh.YData) == numel(histData)
                bh.YData = histData;
            else
                warning("Bar data size mismatch at (%d,%d)", circleIdx, i);
            end
            
        end
    end
    frame = getframe(gcf);
    histPatchVelocitiesVideo(:,:,:,frameIdx) = frame2im(frame);

end
writeGifOnDisc(mat2gray(histPatchVelocitiesVideo), sprintf("histogram_velocities_overlay_%s", name), "ToolBox", ToolBox);
end
close(fi)
end
