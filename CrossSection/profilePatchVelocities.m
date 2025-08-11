function profilePatchVelocities(v_profiles_cell, name, locsLabel, M0_ff_img)

ToolBox = getGlobalToolBox;

params = ToolBox.getParams;
exportVideos = params.exportVideos;

% Check sizes
[rows, cols] = size(locsLabel);
assert(isequal(size(v_profiles_cell), [rows, cols]), 'Size of v_profiles_cell must match locsLabel');
i = 0;
numFrames = 0;

while numFrames <= 0
    i = i + 1;
    numFrames = size(v_profiles_cell{i}, 2);

    if i > size(v_profiles_cell, 1) * size(v_profiles_cell, 2)
        warning("Velocity profiles cells are all empty.")
        break
    end

end

tmp = v_profiles_cell{i};
sizeProfiles = size(tmp{i}, 2) * 2/3;

fi = figure("Visible", "on", 'Color', 'w');
imshow(M0_ff_img, []);
axis image
axis off
fi.Position = [200 200 600 600];

hold on;
title(['Velocity Profiles Overlay - ' name]);

% Parameters for Profiles size
profWidth = 40;
profHeight = 30;

% AVG Plot
lines_cell = cell(rows, cols);

for circleIdx = 1:rows

    for i = 1:cols

        if isempty(locsLabel{circleIdx, i}) || isempty(v_profiles_cell{circleIdx, i})
            continue;
        end

        % Get prof data
        profData = v_profiles_cell{circleIdx, i};

        if ~isequal(size(profData, 2), numFrames)
            warning('Expected v_profiles_cell{%d,%d} to be profile size, numFrames', circleIdx, i);
            continue;
        end

        % Calculate AVG profile
        profile = zeros(size(profData{1}));

        for ff = 1:numFrames
            profile = profile + profData{ff};
        end

        profile = profile / numFrames;

        % Compute axes center location
        pos = locsLabel{circleIdx, i}; % pos = [x, y]

        if isempty(pos) || numel(pos) ~= 2
            continue;
        end

        x = pos(1);
        y = pos(2);
        profile = profile / 50; % Normalize by 50 mm/s
        n = numel(profile);
        x_axis = linspace(-sizeProfiles / 2, sizeProfiles / 2, n);

        % Optional Poiseuille fit
        try

            selected_points = find(profile > max(profile(:)) * 0.3);
            warning('off')
            fitObj = fit(x_axis(selected_points)', double(profile(selected_points)'), 'a*(1 - ((x - b)/c)^2)', ...
                'StartPoint', [1, 0, profWidth / 2]);
            warning('on')
            plotData = fitObj(x_axis);
        catch e
            disp(e)
            plotData = profile;
        end

        % Plot profile
        x_plot = x + x_axis;
        y_data = y - profile * profHeight; % Measured data (true profile)
        y_fit = y - plotData(plotData > 0) * profHeight; % Fitted data

        % Plot directly on image (no text, no axes)
        hold on;
        plot(x_plot(plotData > 0), y_fit, 'r', 'LineWidth', 1.5); % red for fit
        lines_cell{circleIdx, i} = plot(x_plot, y_data, 'k', 'LineWidth', 1); % black for true data

        hold off;
    end

end

% Save figure
exportgraphics(gca, fullfile(ToolBox.path_png, sprintf("%s_velocities_profiles_overlay_%s.png", ToolBox.folder_name, name)));

% Time plot (gif)

if exportVideos
    hold on;
    %imshow(M0_ff_img, []);
    frame = getframe(gca);
    [numX_fig, numY_fig, ~] = size(frame.cdata);
    profilePatchVelocitiesVideo = zeros(numX_fig, numY_fig, 3, numFrames, 'single');

    for frameIdx = 1:numFrames
        %fprintf(" %d ",frameIdx);
        for circleIdx = 1:rows

            for i = 1:cols

                if isempty(locsLabel{circleIdx, i}) || isempty(v_profiles_cell{circleIdx, i})
                    continue;
                end

                % Get prof data
                profData = v_profiles_cell{circleIdx, i};

                if ~isequal(size(profData, 2), numFrames)
                    warning('Expected v_profiles_cell{%d,%d} to be profile size, numFrames', circleIdx, i);
                    continue;
                end

                % Get frame profile
                profile = profData{frameIdx};

                % Compute axes center location
                pos = locsLabel{circleIdx, i}; % pos = [x, y]

                if isempty(pos) || numel(pos) ~= 2
                    continue;
                end

                x = pos(1);
                y = pos(2);
                profile = profile / 50; % Normalize by 50 mm/s
                n = numel(profile);
                x_axis = linspace(-sizeProfiles / 2, sizeProfiles / 2, n);

                % Optional Poiseuille fit
                % try

                %     selected_points = find(profile>max(profile(:))*0.3);
                %     fitObj = fit(x_axis(selected_points)', double(profile(selected_points)'), 'a*(1 - ((x - b)/c)^2)', ...
                %         'StartPoint', [1, 0, profWidth/2]);
                %     plotData = fitObj(x_axis);
                % catch e
                %     disp(e)
                % plotData = profile;

                % Plot profile
                x_plot = x + x_axis;
                y_data = y - profile * profHeight; % Measured data (true profile)
                % y_fit  = y - plotData(plotData>0) * profHeight;    % Fitted data

                % Plot directly on image (no text, no axes)
                hold on;
                set(lines_cell{circleIdx, i}, 'XData', x_plot, 'YData', y_data);

                % plot(x_plot(plotData>0), y_fit, 'r', 'LineWidth', 1.5);    % red for fit
                hold off;
            end

        end

        frame = getframe(gca);
        profilePatchVelocitiesVideo(:, :, :, frameIdx) = frame2im(frame);

    end

    writeGifOnDisc(mat2gray(profilePatchVelocitiesVideo), sprintf("profiles_patch_overlay_%s", name), "ToolBox", ToolBox);

end

close(fi);
end
