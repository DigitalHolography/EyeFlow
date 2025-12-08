function circleImages(M0_ff_img, xy_barycenter, A_cell, Q_cell, v_cell, mask, locsLabel, name)

% Get global ToolBox and parameters
ToolBox = getGlobalToolBox;
params = ToolBox.getParams;
saveFigures = params.saveFigures;
exportVideos = params.exportVideos;

% Validate input name
if ~ismember(name, {'artery', 'vein'})
    error('Name must be either ''artery'' or ''vein''');
end

alphaWom = zeros(size(ToolBox.Cache.WomersleyOut),'single');
for i = 1:size(alphaWom, 1)
    for j = 1:size(alphaWom, 2)
        if isstruct(ToolBox.Cache.WomersleyOut{i, j})
            data = ToolBox.Cache.WomersleyOut{i, j}.segments_metrics.MovingWallFixedNu;
            alphaWom(i, j) = data.alpha_n;
        end
    end
end


path_txt = ToolBox.path_txt;
main_folder = ToolBox.folder_name;
[numX, numY, ~] = size(M0_ff_img);
[numCircles, numBranches] = size(A_cell);

if saveFigures
    path_png = ToolBox.path_png;
    path_eps = ToolBox.path_eps;

    % Extract barycenter coordinates and image dimensions
    x_c = xy_barycenter(1);
    y_c = xy_barycenter(2);
    fontsize = round(numX / 40);

    % Initialize video arrays with consistent size
    videoSize = [numX, numY];
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

    if strcmp(name, 'artery')
        cmap = ToolBox.Cache.cmapArtery;
    else
        cmap = ToolBox.Cache.cmapVein;
    end

    textOptions = {"FontWeight", "bold", ...
        "Color", "white", ...
        "FontSize", fontsize, ...
        "BackgroundColor", "black"};

    % Plot cross-section widths
    parfor cIdx = 1:numCircles

        % Create RGB image with mask overlay
        image_RGB = setcmap(M0_ff_img, maskRadius(:, :, cIdx), cmap) + M0_ff_img .* inv_maskRadius(:, :, cIdx);

        % Create all figures at once
        fig1 = figure("Visible", "off");
        imshow(image_RGB);
        fig2 = figure("Visible", "off");
        imshow(image_RGB);
        fig3 = figure("Visible", "off");
        imshow(image_RGB);
        fig4 = figure("Visible", "off");
        imshow(image_RGB);
        fig5 = figure("Visible", "off");
        imshow(image_RGB);

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
                text(new_x, new_y, sprintf("%d", bIdx), textOptions{:});

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

                % Velocity visualization (if data exists)
                if ~isempty(alphaWom(cIdx, bIdx))
                    figure(fig5);
                    text(new_x, new_y, sprintf('%.1f', alphaWom(cIdx, bIdx)), textOptions{:});
                end

            end

        end



        

        % Capture and resize frame for video
        figure(fig1);
        capturedFrame = frame2im(getframe(gca));
        resizedFrame = rescale(imresize(capturedFrame, [numX, numY]));
        vesselD_video(:, :, :, cIdx) = resizedFrame;

        % Export plot
        exportgraphics(gca, fullfile(path_png, 'vesselSegmentImages', 'lumenDiameter', ...
            sprintf("%s_circle_%d_crossSectionWidth_%s.png", main_folder, cIdx, name)));
        exportgraphics(gca, fullfile(path_eps, 'vesselSegmentImages', 'lumenDiameter', ...
            sprintf("%s_circle_%d_crossSectionWidth_%s.eps", main_folder, cIdx, name)));

        % Capture and resize frame for video
        figure(fig2);
        capturedFrame = frame2im(getframe(gca));
        resizedFrame = rescale(imresize(capturedFrame, [numX, numY]));
        vesselNum_video(:, :, :, cIdx) = resizedFrame;

        % Export plot
        exportgraphics(gca, fullfile(path_png, 'vesselSegmentImages', 'vesselSegmentId', ...
            sprintf("%s_circle_%d_Numerotation%sImage.png", main_folder, cIdx, name)));
        exportgraphics(gca, fullfile(path_eps, 'vesselSegmentImages', 'vesselSegmentId', ...
            sprintf("%s_circle_%d_Numerotation%sImage.eps", main_folder, cIdx, name)));

        % Capture and resize frame for video
        figure(fig3);
        capturedFrame = frame2im(getframe(gca));
        resizedFrame = rescale(imresize(capturedFrame, [numX, numY]));
        Q_video(:, :, :, cIdx) = resizedFrame;

        % Export plot
        exportgraphics(gca, fullfile(path_png, 'vesselSegmentImages', 'bloodVolumeRate', ...
            sprintf("%s_circle_%d_BVR_%s.png", main_folder, cIdx, name)));
        exportgraphics(gca, fullfile(path_eps, 'vesselSegmentImages', 'bloodVolumeRate', ...
            sprintf("%s_circle_%d_BVR_%s.eps", main_folder, cIdx, name)));

        % Capture and resize frame for video
        figure(fig4);
        capturedFrame = frame2im(getframe(gca));
        resizedFrame = rescale(imresize(capturedFrame, [numX, numY]));
        velocity_video(:, :, :, cIdx) = resizedFrame;

        % Export plot
        exportgraphics(gca, fullfile(path_png, 'vesselSegmentImages', 'velocity', ...
            sprintf("%s_circle_%d_velocity%sImage.png", main_folder, cIdx, name)));
        exportgraphics(gca, fullfile(path_eps, 'vesselSegmentImages', 'velocity', ...
            sprintf("%s_circle_%d_velocity%sImage.eps", main_folder, cIdx, name)));

        % Capture and resize frame for video
        figure(fig5);
        capturedFrame = frame2im(getframe(gca));
        resizedFrame = rescale(imresize(capturedFrame, [numX, numY]));
        alpha_video(:, :, :, cIdx) = resizedFrame;

        % Export plot
        exportgraphics(gca, fullfile(path_png, 'vesselSegmentImages', 'alphaWomersley', ...
            sprintf("%s_circle_%d_alpha%sImage.png", main_folder, cIdx, name)));
        exportgraphics(gca, fullfile(path_eps, 'vesselSegmentImages', 'alphaWomersley', ...
            sprintf("%s_circle_%d_alpha%sImage.eps", main_folder, cIdx, name)));


    end

    % Export videos if enabled
    if exportVideos
        writeGifOnDisc(vesselD_video, sprintf('sectionsWidth_%s', name), 0.15, 10);
        writeGifOnDisc(vesselNum_video, sprintf('num_%s', name), 0.15, 10);
        writeGifOnDisc(Q_video, sprintf('BVR_%s', name), 0.15, 10);
        writeGifOnDisc(velocity_video, sprintf('velocity_%s', name), 0.15, 10);
        writeGifOnDisc(alpha_video, sprintf('alpha_womersley_%s', name), 0.15, 10);
    end

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

end
