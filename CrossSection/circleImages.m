function circleImages(M0_ff_img, xy_barycenter, A_cell, Q_cell, v_cell, mask, locsLabel, name)

% Get global toolbox and parameters
ToolBox = getGlobalToolBox;
params = ToolBox.getParams;
exportVideos = params.exportVideos;

% Validate input name
if ~ismember(name, {'Artery', 'Vein'})
    error('Name must be either ''Artery'' or ''Vein''');
end

path_png = ToolBox.path_png;
path_eps = ToolBox.path_eps;
path_txt = ToolBox.path_txt;
main_folder = ToolBox.main_foldername;
initial = name(1);

% Extract barycenter coordinates and image dimensions
x_c = xy_barycenter(1);
y_c = xy_barycenter(2);
[numX, numY, ~] = size(M0_ff_img);
[numCircles, numBranches] = size(A_cell);

% Initialize video arrays with consistent size
videoSize = [465, 465];
vesselD_video = zeros([videoSize, 3, numCircles]);
vesselNum_video = zeros([videoSize, 3, numCircles]);
Q_video = zeros([videoSize, 3, numCircles]);
velocity_video = zeros([videoSize, 3, numCircles]);

% Precompute mask data
maskMat = zeros(numX, numY, numCircles, numBranches);

parfor cIdx = 1:numCircles

    for bIdx = 1:numBranches
        maskMat(:, :, cIdx, bIdx) = mask{cIdx, bIdx};
    end

end

maskRadius = squeeze(sum(maskMat, 4));
inv_maskRadius = ~maskRadius;

% Get appropriate colormap
cmap = ToolBox.(['cmap' name]);
textOptions = {"FontWeight", "bold", ...
    "Color", "white", ...
    "BackgroundColor", "black", ...
    "FontSize", 14};

% Plot cross-section widths
parfor cIdx = 1:numCircles

    % Create RGB image with mask overlay
    image_RGB = setcmap(M0_ff_img, maskRadius(:, :, cIdx), cmap) + M0_ff_img .* inv_maskRadius(:, :, cIdx);

    % Create all figures at once
    fig1 = figure("Visible", "off");
    imagesc(image_RGB);
    fig2 = figure("Visible", "off");
    imagesc(image_RGB);
    fig3 = figure("Visible", "off");
    imagesc(image_RGB);
    fig4 = figure("Visible", "off");
    imagesc(image_RGB);

    for bIdx = 1:numBranches

        if ~isempty(A_cell{cIdx, bIdx}) && ~isnan(A_cell{cIdx, bIdx})

            % Add labels for cross-section widths
            x_l = locsLabel{cIdx, bIdx}(1);
            y_l = locsLabel{cIdx, bIdx}(2);

            c_x = (x_l - x_c) / sqrt((x_l - x_c) ^ 2 + (y_l - y_c) ^ 2);
            c_y = (y_l - y_c) / sqrt((x_l - x_c) ^ 2 + (y_l - y_c) ^ 2);

            new_x = x_l + c_x * numX / 20 - 15;
            new_y = y_l + c_y * numY / 20 - 15;

            % Width visualization
            figure(fig1);
            crossSectionWidth = 2 * sqrt(squeeze(A_cell{cIdx, bIdx}) / pi) * 1000;
            text(new_x, new_y, sprintf("%d µm", round(crossSectionWidth, 0)), textOptions{:});

            % Numerotation visualization
            figure(fig2);
            text(new_x, new_y, sprintf("%s%d", initial, bIdx), textOptions{:});

            % BVR visualization (if data exists)
            if ~isempty(Q_cell{cIdx, bIdx})
                figure(fig3);
                Q = Q_cell{cIdx, bIdx};
                text(new_x, new_y, sprintf('%.1f', mean(Q, 2)), textOptions{:});
            end

            % Velocity visualization (if data exists)
            if ~isempty(v_cell{cIdx, bIdx})
                figure(fig4);
                v = v_cell{cIdx, bIdx};
                text(new_x, new_y, sprintf('%.1f', mean(v, 2)), textOptions{:});
            end

        end

    end

    % Capture and resize frame for video
    figure(fig1);
    capturedFrame = frame2im(getframe(gca));
    resizedFrame = imresize(capturedFrame, [465, 465]);
    vesselD_video(:, :, :, cIdx) = rescale(resizedFrame);

    % Export plot
    exportgraphics(gca, fullfile(path_png, 'crossSectionsAnalysis', 'sectionsImages', 'widths', ...
        sprintf("%s_circle_%d_crossSectionWidth_%s.png", main_folder, cIdx, name)));
    exportgraphics(gca, fullfile(path_eps, 'crossSectionsAnalysis', 'sectionsImages', 'widths', ...
        sprintf("%s_circle_%d_crossSectionWidth_%s.eps", main_folder, cIdx, name)));

    % Capture and resize frame for video
    figure(fig2);
    capturedFrame = frame2im(getframe(gca));
    resizedFrame = imresize(capturedFrame, [465, 465]);
    vesselNum_video(:, :, :, cIdx) = rescale(resizedFrame);

    % Export plot
    exportgraphics(gca, fullfile(path_png, 'crossSectionsAnalysis', 'sectionsImages', 'num', ...
        sprintf("%s_circle_%d_Numerotation%sImage.png", main_folder, cIdx, name)));
    exportgraphics(gca, fullfile(path_eps, 'crossSectionsAnalysis', 'sectionsImages', 'num', ...
        sprintf("%s_circle_%d_Numerotation%sImage.eps", main_folder, cIdx, name)));

    % Capture and resize frame for video
    figure(fig3);
    capturedFrame = frame2im(getframe(gca));
    resizedFrame = imresize(capturedFrame, [465, 465]);
    Q_video(:, :, :, cIdx) = rescale(resizedFrame);

    % Export plot
    exportgraphics(gca, fullfile(path_png, 'crossSectionsAnalysis', 'sectionsImages', 'bvr', ...
        sprintf("%s_circle_%d_BVR_%s.png", main_folder, cIdx, name)));
    exportgraphics(gca, fullfile(path_eps, 'crossSectionsAnalysis', 'sectionsImages', 'bvr', ...
        sprintf("%s_circle_%d_BVR_%s.eps", main_folder, cIdx, name)));

    % Capture and resize frame for video
    figure(fig4);
    capturedFrame = frame2im(getframe(gca));
    resizedFrame = imresize(capturedFrame, [465, 465]);
    velocity_video(:, :, :, cIdx) = rescale(resizedFrame);

    % Export plot
    exportgraphics(gca, fullfile(path_png, 'crossSectionsAnalysis', 'sectionsImages', 'vel', ...
        sprintf("%s_circle_%d_velocity%sImage.png", main_folder, cIdx, name)));
    exportgraphics(gca, fullfile(path_eps, 'crossSectionsAnalysis', 'sectionsImages', 'vel', ...
        sprintf("%s_circle_%d_velocity%sImage.eps", main_folder, cIdx, name)));
end

% Create and save table of vessel widths
widthtable = table;
widthtable.Label = (1:numBranches)';

for cIdx = 1:numCircles
    crossSectionWidth = zeros(numBranches, 1);

    for bIdx = 1:numBranches

        if ~isempty(A_cell{cIdx, bIdx})
            crossSectionWidth(bIdx) = 2 * sqrt(squeeze(A_cell{cIdx, bIdx}) / pi) * 1000;
        end

    end

    widthtable.(sprintf("circle_%d", cIdx)) = crossSectionWidth;
end

writetable(widthtable, fullfile(path_txt, strcat(main_folder, '_Vessels_Widths_Table.txt')));

% Export videos if enabled
if exportVideos
    writeGifOnDisc(vesselD_video, sprintf('sectionsWidth_%s', name), 0.15, 10);
    writeGifOnDisc(vesselNum_video, sprintf('num_%s', name), 0.15, 10);
    writeGifOnDisc(Q_video, sprintf('BVR_%s', name), 0.15, 10);
    writeGifOnDisc(velocity_video, sprintf('velocity_%s', name), 0.15, 10);
end

end
