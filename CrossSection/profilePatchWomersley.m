function profilePatchWomersley(v_profiles_cell, name, locsLabel, M0_ff_img)

ToolBox = getGlobalToolBox;

% Check sizes and extract numFrames from first non empty profile data in input
[rows, cols] = size(locsLabel);
assert(isequal(size(v_profiles_cell), [rows, cols]), 'Size of v_profiles_cell must match locsLabel');
ind = 0;
numFrames = 0;

while numFrames <= 0
    ind = ind + 1;
    numFrames = size(v_profiles_cell{ind}, 2);

    if ind > size(v_profiles_cell, 1) * size(v_profiles_cell, 2)
        warning("Velocity profiles cells are all empty.")
        break
    end

end

% Extract cardiac frequency and corresponding indices with a margin
cardiac_frequency = ToolBox.Cache.HeartBeatFFT; % in Hz

f = fft_freq_vector(ToolBox.fs * 1000 / ToolBox.stride, numFrames);

[~, cardiac_idx] = min(abs(f - cardiac_frequency));

margin_ = round(0.1 * (cardiac_idx - numFrames / 2)); % +- 10 % of Heartrate
cardiac_idxs = cardiac_idx + (-margin_:margin_);
cardiac_idxs(cardiac_idxs > numFrames) = [];
cardiac_idxs(cardiac_idxs < 1) = [];

% Start plotting
tmp = v_profiles_cell{ind};
sizeProfiles = size(tmp{ind}, 2) * 2/3;

fi = figure("Visible", "off", 'Color', 'w');
imshow(M0_ff_img, []);
axis image
axis off
fi.Position = [200 200 600 600];

hold on;
title(['Womersley Profiles Overlay - ' name]);

% Parameters for Profiles size
% profWidth = 40;
profHeight = 30;

% AVG Plot
% lines_cell = cell(rows, cols);

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

        % Calculate FFT of the time dependent profile
        profile_time = zeros(length(profData{1}), numFrames);

        for ff = 1:numFrames
            profile_time(:, ff) = profData{ff};
        end

        profile_ft = fftshift(fft(profile_time, [], 2), 2);

        % Calculate the complex Wom profile to plot

        profile_Wom = mean(profile_ft(:, cardiac_idxs), 2);

        profile_Wom = profile_Wom / mean(profile_Wom);

        % Compute axes center location
        pos = locsLabel{circleIdx, i}; % pos = [x, y]

        if isempty(pos) || numel(pos) ~= 2
            continue;
        end

        x = pos(1);
        y = pos(2);
        profile_Wom = profile_Wom / 5; % Normalize by 50
        n = numel(profile_Wom);
        x_axis = linspace(-sizeProfiles / 2, sizeProfiles / 2, n);

        % Plot profile
        x_plot = x + x_axis;
        y_data_r = y - real(profile_Wom) * profHeight; % Measured data (true profile)
        y_data_i = y - imag(profile_Wom) * profHeight; % Measured data (true profile)

        % Plot directly on image (no text, no axes)
        hold on;
        plot(x_plot, y_data_r, 'b', 'LineWidth', 1); % blue for real
        plot(x_plot, y_data_i, 'r', 'LineWidth', 1); % red for imag

        hold off;
    end

end

% Save figure
exportgraphics(gca, fullfile(ToolBox.path_png, sprintf("%s_velocities_womersley_profiles_overlay_%s.png", ToolBox.folder_name, name)));

close(fi);
end
