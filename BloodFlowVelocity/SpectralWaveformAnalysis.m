function [fft_c, fundamental, valid_harmonics, f] = SpectralWaveformAnalysis(signal, numSys, name)
% Spectral Analysis
% Perform spectral analysis on the original signal
% Zero-pad the signal for better frequency resolution

ToolBox = getGlobalToolBox;
params = ToolBox.getParams;
saveFigures = params.saveFigures;
numFrames = length(signal);
fs = 1 / (ToolBox.stride / ToolBox.fs / 1000); % Sampling frequency in Hz
duration = numFrames * ToolBox.stride / ToolBox.fs / 1000;
estimated_fundamental = numSys / duration;

if strcmp(name, "v_artery")
    cDark = [1 0 0]; % Dark color for peaks
    fig_name = 'ArterialSpectralAnalysis';
elseif strcmp(name, "v_vein")
    cDark = [0 0 1]; % Dark color for peaks
    fig_name = 'VenousSpectralAnalysis';
else
    cDark = [0 0 0]; % Default color
    fig_name = 'SpectralAnalysis';
end

% --- APPLY HAMMING WINDOW TO TIME-DOMAIN SIGNAL ---
hamming_win = hamming(numFrames)'; % Create a Hamming window of the same length as the original signal
windowed_signal = signal .* hamming_win; % Apply the window (element-wise multiplication)

% Calculate coherent gain of the Hamming window for amplitude compensation
coherent_gain = mean(hamming_win);

% Zero-pad the WINDOWED signal
N = 2; % Padding factor
padded_signal = padarray(windowed_signal, [0, numFrames * N]);

% Frequency vector (show only positive frequencies since signal is real)
f = fft_freq_vector(fs, length(padded_signal), true);
fft_c = fft(padded_signal);

% --- COMPENSATE FOR WINDOW AMPLITUDE LOSS ---
fft_mag = abs(fft_c) / coherent_gain; % Divide by coherent gain

fft_mag = fft_mag(1:length(f)); % Take only positive frequencies
fft_mag = fft_mag / max(fft_mag); % Normalize to [0,1]

% Improved peak detection with minimum prominence threshold
fundamental_peak = findpeaks(fft_mag, ...
    'MinPeakDistance', estimated_fundamental * 0.6, ...
    'SortStr', 'descend', 'NPeaks', 1);
min_prominence = fundamental_peak(1) * 0.1; % 10 % of maximum magnitude as minimum prominence

% Improved peak detection with minimum prominence threshold
[s_peaks, s_idx] = findpeaks(fft_mag, ...
    'MinPeakDistance', estimated_fundamental * 0.8, ...
    'MinPeakProminence', min_prominence, ...
    'SortStr', 'descend');
s_locs = f(s_idx);
numFreq = length(s_locs);

% Sort peaks by descending magnitude
[s_locs, idx] = sort(s_locs, 'ascend');
s_peaks = s_peaks(idx);
s_idx = s_idx(idx);
m_harmonics = length(s_locs);

% Calculate and plot harmonic frequencies if fundamental is detected
if numFreq >= 2
    valid_harmonics = []; % Initialize valid harmonics array
    fundamental = s_locs(1);
    harmonics = fundamental * (2:length(s_peaks)); % Up to m th harmonic

    % For each harmonic, find the closest local maximum (peak) in the spectrum
    valid_harmonics(1) = fundamental;
    [~, locs] = findpeaks(fft_mag);

    for h = harmonics
        % Define a small window around the harmonic frequency (e.g., ±0.05*fundamental)
        window = fundamental * 0.48;
        idx_window = find(f >= h - window & f <= h + window);

        if isempty(idx_window)
            continue;
        end

        local_locs = locs(locs >= idx_window(1) & locs <= idx_window(end));

        if ~isempty(local_locs)
            % Choose the peak closest to the harmonic frequency
            [~, min_idx] = min(abs(f((local_locs)) - h));
            valid_harmonics(end + 1) = f((local_locs(min_idx))); %#ok<AGROW>
        end

    end

    valid_harmonics = valid_harmonics(valid_harmonics <= fs / 2);

elseif numFreq == 1
    fundamental = s_locs(1);
    valid_harmonics = fundamental;
    valid_harmonics = valid_harmonics(valid_harmonics <= fs / 2);

else
    fundamental = NaN;
    valid_harmonics = [];
end

% Save to ToolBox
ToolBox.Cache.harmonics = valid_harmonics;

% Add heart rate information if fundamental is in typical range
if ~isempty(s_locs)
    freqs = round(s_locs / fundamental);
    hr = (freqs' \ s_locs');
    freq_estim = freqs * hr;
    residus = s_locs - freq_estim;
    RMSE = sqrt(mean(residus .^ 2));
    hr_se = RMSE / sqrt(numFreq);

    ToolBox.Output.add('HeartBeatFFT', 60 * hr, 'bpm', 60 * hr_se, 'h5path', '/Artery/HeartBeatFFT');
    ToolBox.Cache.HeartBeatFFT = hr; % Save heart rate to cache in Hz
    ToolBox.Cache.HeartBeatFFTSTE = hr_se; % Save heart rate standard error to cache in Hz
else
    ToolBox.Output.add('HeartBeatFFT', NaN, 'bpm', NaN, 'h5path', '/Artery/HeartBeatFFT');
    ToolBox.Cache.HeartBeatFFT = NaN; % Save heart rate to cache in Hz
    ToolBox.Cache.HeartBeatFFTSTE = NaN; % Save heart rate standard error to cache in Hz
end

if saveFigures
    % --- MAGNITUDE SPECTRUM ANALYSIS ---
    % Create figure for spectral analysis
    hFig = figure('Visible', 'off', 'Color', 'w');

    % Main plot with improved styling
    plot(f, fft_mag, 'k', 'LineWidth', 2);
    hold on;
    grid on;

    % Highlight detected peaks with annotations
    if ~isempty(s_peaks)
        scatter(s_locs, s_peaks, 100, 'filled', 'MarkerFaceColor', cDark, 'MarkerEdgeColor', 'k');
    end

    % Configure axes
    axis_height = 1.3 * fundamental_peak; % Set y-axis limit slightly above the highest peak
    axis tight;
    axT = axis;
    axis([axT(1), 10, 0, axis_height]);
    pbaspect([1.618 1 1]);
    set(gca, 'LineWidth', 2, 'FontSize', 12);

    % Annotate the top peaks
    for k = 1:numFreq
        % Plot fundamental location
        xline(s_locs(k), '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1, ...
            'Label', sprintf('%.0f×', s_locs(k) / fundamental), ...
            'LabelOrientation', 'horizontal', ...
            'LabelVerticalAlignment', 'bottom');

        text(s_locs(k), s_peaks(k) + 0.1 * axis_height, ...
            sprintf('%.2f', s_locs(k)), ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'bottom', ...
            'FontSize', 10, ...
            'BackgroundColor', 'w', ...
            'EdgeColor', 'none');
    end

    % Labels and title
    xlabel('Frequency (Hz)', 'FontSize', 14);
    ylabel('Normalized Magnitude', 'FontSize', 14);

    if ~isempty(s_locs)

        annotation('textbox', [0.5 0.6 0.2 0.1], ...
            'String', sprintf('HR : %.1f BPM ± %.1f', 60 * hr, 60 * hr_se), ... % Convert Hz to BPM
            'FitBoxToText', 'on', ...
            'BackgroundColor', 'w', ...
            'EdgeColor', 'none', ...
            'FontSize', 12);
    end

    % Save Results
    exportgraphics(hFig, fullfile(ToolBox.path_png, ...
        sprintf("%s_%s_%s.png", ToolBox.folder_name, fig_name, name)), ...
        'Resolution', 300);

    % Close the figure if not needed
    if ~strcmpi(get(hFig, 'Visible'), 'on')
        close(hFig);
    end

    % --- PHASE SPECTRUM ANALYSIS ---

    % Main plot with improved styling
    fft_angle = angle(fft_c);
    fft_angle = fft_angle(1:length(f)); % Take only positive frequencies

    % Create figure for spectral analysis
    hFig_angle = figure('Visible', 'off', 'Color', 'w');
    plot(f, fft_angle, 'k', 'LineWidth', 2);
    hold on;
    grid on;

    % Highlight detected peaks with annotations
    if ~isempty(s_peaks)
        scatter(s_locs, fft_angle(s_idx), 100, 'filled', 'MarkerFaceColor', cDark, 'MarkerEdgeColor', 'k');

        % Annotate the top peaks
        for k = 1:length(s_peaks)
            text(s_locs(k), fft_angle(s_idx(k)) + 0.3, ...
                sprintf('%.2f', s_locs(k)), ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'bottom', ...
                'FontSize', 10, ...
                'BackgroundColor', 'w', ...
                'EdgeColor', 'none');
        end

    end

    % Configure axes
    axis([0, 10, -pi, pi]);

    xlabel('Frequency (Hz)', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('Phase (rad)', 'FontSize', 14, 'FontWeight', 'bold');
    pbaspect([1.618 1 1]);
    set(gca, 'LineWidth', 1.5, 'FontSize', 12);

    % Save Results
    exportgraphics(hFig_angle, fullfile(ToolBox.path_png, ...
        sprintf("%s_%s_Phase_%s.png", ToolBox.folder_name, fig_name, name)), ...
        'Resolution', 300);

    % Close the figure if not needed
    if ~strcmpi(get(hFig_angle, 'Visible'), 'on')
        close(hFig_angle);
    end

end

end
