function [fft_c, fundamental, valid_harmonics] = SpectralWaveformAnalysis(signal, m_harmonics, name)
% Spectral Analysis
% Perform spectral analysis on the original signal
% Zero-pad the signal for better frequency resolution

ToolBox = getGlobalToolBox;
numFrames = length(signal);
fs = 1 / (ToolBox.stride / ToolBox.fs / 1000);
t = linspace(0, numFrames / fs, numFrames);

cDark = [1 0 0];
cLight = [1 0.5 0.5];

% Set parameters based on signal type
if strcmpi(name, "bvr")
    y_label = 'Flow Rate (µL/min)';
    unit = 'µL/min';
    isBVR = true;
else
    y_label = 'Velocity (mm/s)';
    unit = 'mm/s';
    isBVR = false;
end

N = 16; % Padding factor
padded_signal = padarray(signal, [0 numFrames * N]); % Zero-padding for interpolation in frequency domain

% Frequency vector (show only positive frequencies since signal is real)
f = linspace(0, fs / 2, (N * numFrames) + 1);
fft_c = fft(padded_signal);
fft_mag = abs(fft_c);
fft_mag = fft_mag(1:length(f)); % Take only positive frequencies
fft_mag = fft_mag / max(fft_mag); % Normalize to [0,1]

% Improved peak detection with minimum prominence threshold
min_prominence = 0.1; % 10 % of maximum magnitude as minimum prominence
f_indx = f>0.5; % search only in f > 0.5 Hz
[s_peaks, s_locs] = findpeaks(fft_mag(f_indx), f(f_indx), ...
    'MinPeakProminence', min_prominence, ...
    'SortStr', 'descend', ...
    'NPeaks', 5); % Find up to 5 most significant peaks

% Sort peaks by descending magnitude
[s_peaks, idx] = sort(s_peaks, 'descend');
s_locs = s_locs(idx);

% Create figure for spectral analysis
hFig = figure('Visible', 'on', 'Color', 'w');

% Main plot with improved styling
plot(f, fft_mag, 'k', 'LineWidth', 2);
hold on;
grid on;

% Highlight detected peaks with annotations
if ~isempty(s_peaks)
    scatter(s_locs, s_peaks, 100, 'filled', 'MarkerFaceColor', cDark, 'MarkerEdgeColor', 'k');

    % Annotate the top 3 peaks
    for k = 1:min(5, length(s_peaks))
        text(s_locs(k), s_peaks(k) * 1.2, ...
            sprintf('%.2f Hz', s_locs(k)), ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'bottom', ...
            'FontSize', 8, ...
            'BackgroundColor', 'w', ...
            'EdgeColor', 'k');
    end

end

% Calculate and plot harmonic frequencies if fundamental is detected
if length(s_locs) >= 1
    fundamental = s_locs(1);
    harmonics = fundamental * (2:m_harmonics); % Up to m th harmonic
    % For each harmonic, find the closest local maximum (peak) in the spectrum
    valid_harmonics = [];
    [pks, locs] = findpeaks(fft_mag,'MinPeakProminence',min_prominence/4,'NPeaks', 15);
    for h = harmonics
        % Define a small window around the harmonic frequency (e.g., ±0.05*fundamental)
        window = fundamental * 1.18;
        idx_window = find(f >= h - window & f <= h + window);
        if isempty(idx_window)
            continue;
        end
        local_locs = locs(locs >= idx_window(1) & locs <=idx_window(end));
        if ~isempty(local_locs)
            % Choose the peak closest to the harmonic frequency
            [~, min_idx] = min(abs(f((local_locs)) - h));
            valid_harmonics(end+1) = f((local_locs(min_idx))); %#ok<AGROW>
        end
    end
    valid_harmonics = valid_harmonics(valid_harmonics <= fs / 2);

    % Plot harmonic locations
    for h = valid_harmonics
        xline(h, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1, ...
            'Label', sprintf('%.0f×', h / fundamental), ...
            'LabelOrientation', 'horizontal', ...
            'LabelVerticalAlignment', 'bottom');
    end

end

% Save to ToolBox

ToolBox.Cache.list.harmonics = [fundamental valid_harmonics];

% Configure axes
axis tight;
axT = axis;
axis padded;
axP = axis;
axis([axT(1), 10, axP(3) - 0.1 * (axP(4) - axP(3)), axP(4) + 0.1 * (axP(4) - axP(3))]);

xlabel('Frequency (Hz)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Normalized Magnitude', 'FontSize', 14, 'FontWeight', 'bold');
pbaspect([1.618 1 1]);
set(gca, 'LineWidth', 1.5, 'FontSize', 12);

% Add heart rate information if fundamental is in typical range (0.5-3 Hz)
if ~isempty(s_locs) && s_locs(1) >= 0.5 && s_locs(1) <= 3
    hr = s_locs(1) * 60; % Convert Hz to BPM
    annotation('textbox', [0.5 0.6 0.2 0.1], ...
        'String', sprintf('Estimated HR: %.1f BPM', hr), ...
        'FitBoxToText', 'on', ...
        'BackgroundColor', 'w', ...
        'EdgeColor', 'k', ...
        'FontSize', 12);
end

% Save Results
exportgraphics(hFig, fullfile(ToolBox.path_png, ...
    sprintf("%s_ArterialSpectralAnalysis_%s.png", ToolBox.folder_name, name)), ...
    'Resolution', 300);

end
