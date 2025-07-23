function histogramPatchVelocities(histo_v_cell, name, locsLabel, maskLabel, M0_ff_img)

ToolBox = getGlobalToolBox;

params = ToolBox.getParams;
exportVideos = params.exportVideos;

% Check sizes
[rows, cols] = size(locsLabel);
assert(isequal(size(histo_v_cell), [rows, cols]), 'Size of histo_v_cell must match locsLabel');

numFrames = size(histo_v_cell{1},2);
fi = figure();
imshow(M0_ff_img, []);
hold on;
title(['Velocity Histograms Overlay - ' name]);

% Parameters for histogram size
histWidth = 40;
histHeight = 30;
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
        pos = locsLabel{circleIdx, i};  % pos = [x, y]
        if isempty(pos) || numel(pos) ~= 2
            continue;
        end
        x = pos(1);
        y = pos(2);

        % Create new axes for histogram at position
        ax = axes('Position', ...
            [(x - histWidth/2) / size(M0_ff_img, 2), ...
             1 - (y + histHeight/2) / size(M0_ff_img, 1), ...
             histWidth / size(M0_ff_img, 2), ...
             histHeight / size(M0_ff_img, 1)]);

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

outputDir = fullfile(ToolBox.path_png, 'local');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end
% Save figure
saveas(gcf, fullfile(outputDir, ...
    sprintf("%s_histogram_velocities_overlay_%s.png", ToolBox.folder_name, name)));

%% (GIF)
if exportVideos
hold on;

histPatchVelocitiesVideo = zeros(1124,1255, 3, numFrames,'single');
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

end
