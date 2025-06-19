function ArterialWaveformAnalysis(signal, systolesIndexes, numInterp, name)
% ARTERIALWAVEFORMANALYSIS Analyzes arterial waveform signals (velocity or blood volume rate)
%
% Inputs:
%   signal - Input waveform signal
%   systolesIndexes - Indices of systolic peaks
%   numInterp - Number of interpolation points
%   name - Signal type ('bvr' for blood volume rate, otherwise velocity)

%% Initial Setup
ToolBox = getGlobalToolBox;
numFrames = length(signal);
fs = 1 / (ToolBox.stride / ToolBox.fs / 1000);
t = linspace(0, numFrames / fs, numFrames);

cDark = [0.75 0 0];
cLight = [1 0.5 0.5];

% Set parameters based on signal type
if strcmpi(name, "bvr")
    y_label = 'Blood Volume Rate (µL/min)';
    folder = 'crossSectionsAnalysis';
    unit = 'µL/min';
    isBVR = true;
else
    y_label = 'Velocity (mm/s)';
    folder = 'bloodFlowVelocity';
    unit = 'mm/s';
    isBVR = false;
end

%% Signal Preprocessing
try
    % Apply wavelet denoising if possible
    signal = double(wdenoise(signal, 4, ...
        'Wavelet', 'sym4', ...
        'DenoisingMethod', 'Bayes', ...
        'ThresholdRule', 'Median'));
catch
    signal = double(signal);
end

%% Cycle Analysis
[one_cycle_signal, avgLength] = interpSignal(signal, systolesIndexes, numInterp);

% Create time vector for one cycle
dt = (t(2) - t(1));
pulseTime = linspace(0, dt * avgLength, numInterp);

%% Feature Detection
% Adaptive peak detection parameters
min_peak_height = max(one_cycle_signal) * 0.3; % 30 % of max as threshold
min_peak_distance = floor(length(one_cycle_signal) / 4); % 1/4 cycle minimum

% Find peaks with improved parameters
[peaks, locs_peaks] = findpeaks(one_cycle_signal, ...
    'MinPeakHeight', min_peak_height, ...
    'MinPeakDistance', min_peak_distance, ...
    'NPeaks', 2, ... % Limit to 2 peaks max
    'SortStr', 'descend'); % Sort by amplitude

% Initialize all output variables
systoleDuration = NaN;
diastoleDuration = NaN;
systolicUpstroke = peaks(1) - one_cycle_signal(1);
systolicDownstroke = NaN;
diastolicRunoff = NaN;
notch = NaN;
locs_notch = NaN;

% Detect global minimum at end of signal (last 25% of cycle)
endSegment = floor(0.75 * numInterp):numInterp;
[endMinVal, endMinLoc] = min(one_cycle_signal(endSegment));
endMinLoc = endMinLoc + endSegment(1) - 1; % Adjust index

% Detect dicrotic notch if two peaks found
if length(peaks) > 1
    [notch, locs_notch] = min(one_cycle_signal(locs_peaks(1):locs_peaks(2)));
    locs_notch = locs_notch + locs_peaks(1) - 1;

    % Only consider valid notch (significant difference from diastolic peak)
    if (peaks(2) - notch) > peaks(1) * 0.1 % 10 % threshold
        systolicDownstroke = peaks(1) - notch;
        diastolicRunoff = notch - one_cycle_signal(end); % End of cycle

        % Calculate durations
        systoleDuration = pulseTime(locs_notch) - pulseTime(1);
        diastoleDuration = pulseTime(end) - pulseTime(locs_notch);
    else
        notch = NaN; % Invalid notch
    end

end

%% Visualization
hFig = figure('Visible', 'on', 'Color', 'w');
hold on;

% Add reference lines and annotations
T_peak = pulseTime(locs_peaks(1));
xline(T_peak, 'k--', sprintf("%.2f s", T_peak), ...
    'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom', 'Color', [0.25 0.25 0.25], 'FontSize', 10);

xline(dt * avgLength, 'k--', sprintf("%.2f s", dt * avgLength), ...
    'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom', 'Color', [0.25 0.25 0.25], 'FontSize', 10);

yline(peaks(1), 'k--', sprintf("%.1f %s", peaks(1), unit), ...
    'LineWidth', 1.5, 'Color', [0.25 0.25 0.25], 'FontSize', 10);

T_Min = pulseTime(endMinLoc);
xline(T_Min, 'k--', sprintf("%.2f s", T_Min), ...
    'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom', 'Color', [0.25 0.25 0.25], 'FontSize', 10);

yline(endMinVal, 'b--', sprintf("%.1f %s", endMinVal, unit), ...
    'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom', 'Color', [0.25 0.25 0.25], 'FontSize', 10);

if ~isnan(notch)
    T_notch = pulseTime(locs_notch);
    xline(T_notch, 'k--', sprintf("%.2f s", T_notch), ...
        'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom', 'Color', [0.25 0.25 0.25], 'FontSize', 10);

    yline(notch, 'k--', sprintf("%.1f %s", notch, unit), ...
        'LineWidth', 1.5, 'Color', [0.25 0.25 0.25], 'FontSize', 10);

    if length(locs_peaks) > 1
        T_diastolic = pulseTime(locs_peaks(2));
        xline(T_diastolic, 'k--', sprintf("%.2f s", T_diastolic), ...
            'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom', 'Color', [0.25 0.25 0.25], 'FontSize', 10);

        yline(peaks(2), 'k--', sprintf("%.1f %s", peaks(2), unit), ...
            'LineWidth', 1.5, 'Color', [0.25 0.25 0.25], 'FontSize', 10);
    end

end

padded_signal = [one_cycle_signal(floor(numInterp / 2) + 1:end), ...
                     one_cycle_signal, ...
                     one_cycle_signal(1:floor(numInterp / 2))];
padded_time = [linspace(- dt * avgLength / 2, -dt, round(numInterp / 2)), ...
                   pulseTime, ...
                   linspace(dt * (avgLength + 1), dt * (avgLength * 3/2), round(numInterp / 2))];
padded_gradient = [gradient(one_cycle_signal(floor(numInterp / 2) + 1:end)), ...
                       gradient(one_cycle_signal), ...
                       gradient(one_cycle_signal(1:floor(numInterp / 2)))];

% Main signal and gradient plots
plot(padded_time, padded_gradient, 'Color', [0.85 0.85 0.85], 'LineWidth', 2);
plot(padded_time, padded_signal, 'Color', [0.85 0.85 0.85], 'LineWidth', 2);
plot(pulseTime, gradient(one_cycle_signal), 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
plot(pulseTime, one_cycle_signal, 'k', 'LineWidth', 2);

% Plot detected features
scatter(pulseTime(locs_peaks), peaks, 100, 'filled', 'MarkerFaceColor', cDark, 'MarkerEdgeColor', 'k');
scatter(pulseTime(endMinLoc), endMinVal, 100, 'filled', 'MarkerFaceColor', cLight, 'MarkerEdgeColor', 'k');

if ~isnan(notch)
    scatter(pulseTime(locs_notch), notch, 100, 'filled', 'MarkerFaceColor', cLight, 'MarkerEdgeColor', 'k');
end

% Configure axes
axis tight;
axT = axis;
axis padded;
axP = axis;
axis([axT(1), axT(2), axP(3) - 0.1 * (axP(4) - axP(3)), axP(4) + 0.1 * (axP(4) - axP(3))]);

ylabel(y_label);
xlabel('Time (s)');
pbaspect([1.618 1 1]);
set(gca, 'LineWidth', 2, 'Box', 'on');

%% Save Results
exportgraphics(hFig, fullfile(ToolBox.path_png, folder, ...
    sprintf("%s_ArterialWaveformAnalysis_%s.png", ToolBox.folder_name, name)), ...
    'Resolution', 300);

%% Spectral Analysis
% Perform spectral analysis on the original signal
% Zero-pad the signal for better frequency resolution
N = 16; % Padding factor
padded_signal = padarray(signal, [0 numFrames * N]); % Zero-padding for interpolation in frequency domain

% Frequency vector (show only positive frequencies since signal is real)
f = linspace(0, fs / 2, (N * numFrames) + 1);
fft_mag = abs(fft(padded_signal));
fft_mag = fft_mag(1:length(f)); % Take only positive frequencies
fft_mag = fft_mag / max(fft_mag); % Normalize to [0,1]

% Improved peak detection with minimum prominence threshold
min_prominence = 0.1; % 10 % of maximum magnitude as minimum prominence
[s_peaks, s_locs] = findpeaks(fft_mag, f, ...
    'MinPeakProminence', min_prominence, ...
    'SortStr', 'descend', ...
    'NPeaks', 5); % Find up to 5 most significant peaks

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
    harmonics = fundamental * (2:5); % Up to 5th harmonic
    valid_harmonics = harmonics(harmonics <= fs / 2);

    % Plot harmonic locations
    for h = valid_harmonics
        xline(h, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1, ...
            'Label', sprintf('%.0f×', h / fundamental), ...
            'LabelOrientation', 'horizontal', ...
            'LabelVerticalAlignment', 'bottom');
    end

end

% Enhanced axis formatting
xlim([0 10]);
ylim([1e-3 1.1]); % Adjusted for normalized magnitudes
xlabel('Frequency (Hz)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Normalized Magnitude (log scale)', 'FontSize', 14, 'FontWeight', 'bold');
title('Power Spectrum with Peak Detection', 'FontSize', 16);
pbaspect([1.618 1 1]);
set(gca, 'LineWidth', 1.5, 'FontSize', 12);

% Add heart rate information if fundamental is in typical range (0.5-3 Hz)
if ~isempty(s_locs) && s_locs(1) >= 0.5 && s_locs(1) <= 3
    hr = s_locs(1) * 60; % Convert Hz to BPM
    annotation('textbox', [0.7 0.8 0.2 0.1], ...
        'String', sprintf('Estimated HR: %.1f BPM', hr), ...
        'FitBoxToText', 'on', ...
        'BackgroundColor', 'w', ...
        'EdgeColor', 'k', ...
        'FontSize', 12);
end

%% Save Results
exportgraphics(hFig, fullfile(ToolBox.path_png, folder, ...
    sprintf("%s_ArterialSpectralAnalysis_%s.png", ToolBox.folder_name, name)), ...
    'Resolution', 300);

%% Export to JSON (only for velocity signals)
if ~isBVR
    ToolBox.Outputs.add('SystoleDuration', systoleDuration, 's');
    ToolBox.Outputs.add('DiastoleDuration', diastoleDuration, 's');
    ToolBox.Outputs.add('SystolicUpstroke', systolicUpstroke, unit);
    ToolBox.Outputs.add('SystolicDownstroke', systolicDownstroke, unit);
    ToolBox.Outputs.add('DiastolicRunoff', diastolicRunoff, unit);
end

% Close the figure if not needed
if ~strcmpi(get(hFig, 'Visible'), 'on')
    close(hFig);
end

end
