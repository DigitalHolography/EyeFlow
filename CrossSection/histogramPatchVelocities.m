function histogramPatchVelocities(histo_v_cell, name, locsLabel, maskLabel, M0_ff_img)

ToolBox = getGlobalToolBox;
% Check sizes
[rows, cols] = size(locsLabel);
assert(isequal(size(histo_v_cell), [rows, cols]), 'Size of histo_v_cell must match locsLabel');

figure(45);
imshow(M0_ff_img, []);
hold on;
title(['Velocity Histograms Overlay - ' name]);

% Parameters for histogram size
histWidth = 40;
histHeight = 30;

for circleIdx = 1:rows
    for i = 1:cols
        if isempty(locsLabel{circleIdx, i}) || isempty(histo_v_cell{circleIdx, i})
            continue;
        end

        % Get histogram data
        histData = histo_v_cell{circleIdx, i};
        if numel(histData) ~= 2
            warning('Expected histo_v_cell{%d,%d} to be {counts, edges}', circleIdx, i);
            continue;
        end

        counts = histData{1};
        edges = histData{2};
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
        bar(ax, edges(1:end-1), counts, 'histc');
        ax.XTick = [];
        ax.YTick = [];
        ax.Box = 'off';
        ax.Color = 'none';
        ax.XColor = 'none';
        ax.YColor = 'none';
    end
end


outputDir = fullfile(ToolBox.path_png, 'crossSectionsAnalysis');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end
% Save figure
saveas(gcf, fullfile(outputDir, ...
    sprintf("%s_velocities_histogram_overlay_%s.png", ToolBox.folder_name, name)));
close(45);
end
