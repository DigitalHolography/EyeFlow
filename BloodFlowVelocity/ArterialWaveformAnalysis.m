function one_cycle_signal = ArterialWaveformAnalysis(signal, systolesIndexes, numInterp, name)
% ARTERIALWAVEFORMANALYSIS Analyzes arterial waveform signals (velocity or Flow Rate)
%
% Inputs:
%   signal - Input waveform signal
%   systolesIndexes - Indices of systolic peaks
%   numInterp - Number of interpolation points
%   name - Signal type ('bloodVolumeRate' for Flow Rate, otherwise velocity)

% Initial Setup
ToolBox = getGlobalToolBox;
numFrames = length(signal);
fs = 1 / (ToolBox.stride / ToolBox.fs / 1000);
t = linspace(0, numFrames / fs, numFrames);
numSystoles = length(systolesIndexes);

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

% Spectral Analysis
numHarmonics = 6;
[fft_c, ~, valid_harmonics, ~] = SpectralWaveformAnalysis(signal, numSystoles, numHarmonics, name);

% % Signal Preprocessing
% try
%     % Apply wavelet denoising if possible
%     signal = double(wdenoise(signal, 4, ...
%         'Wavelet', 'sym4', ...
%         'DenoisingMethod', 'Bayes', ...
%         'ThresholdRule', 'Median'));
% catch
%     signal = double(signal);
% end

zeroPadLength = length(fft_c);
harmonics_index = round(valid_harmonics / (fs / 2) * zeroPadLength);
fft_abs = abs(fft_c);
fft_angle = angle(fft_c);

% Apply bandpass filter (0.5-15 Hz) as suggested
[b, a] = butter(4, 15 / (fs / 2), 'low');
filtered_signal = filtfilt(b, a, signal);

% Cycle Analysis
[one_cycle_signal, avgLength] = interpSignal(filtered_signal, systolesIndexes, numInterp);
L = length(one_cycle_signal);

% Create time vector for one cycle
dt = (t(2) - t(1));
pulseTime = linspace(0, dt * avgLength, numInterp);

% Create harmonics and combination signal
A_h = fft_abs(harmonics_index);
phi_h = fft_angle(harmonics_index);
composite_signal = A_h .* cos(2 * pi * valid_harmonics .* pulseTime' + phi_h);
composite_signal = composite_signal / max(composite_signal(1, :), [], 'all') * max(signal, [], 'all'); % rescale at initial scale

% Feature Detection

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
dicroticNotchTime = NaN;
diastoleDuration = NaN;
systolicUpstroke = peaks(1) - one_cycle_signal(1);
systolicDownstroke = NaN;
diastolicRunoff = NaN;
% dicroticNotchIndex = NaN;
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
    if (locs_peaks(2) - locs_notch) > L * 0.05 % 5 % threshold
        ToolBox.Output.add("DicroticNotchVisibility", 1, '');
        systolicDownstroke = peaks(1) - notch;
        diastolicRunoff = notch - one_cycle_signal(end); % End of cycle

        % Calculate durations
        dicroticNotchTime = pulseTime(locs_notch) - pulseTime(1);
        % dicroticNotchIndex = notch ./ peaks(1);
        diastoleDuration = pulseTime(end) - pulseTime(locs_notch);
    else
        notch = NaN; % Invalid notch
        ToolBox.Output.add("DicroticNotchVisibility", 0, '');
    end

else
    ToolBox.Output.add("DicroticNotchVisibility", 0, '');
end

% Visualization
hFig = figure('Visible', 'off', 'Color', 'w');
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
        'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom', 'Color', [0.25 0.25 0.25], 'FontSize', 10);

    if length(locs_peaks) > 1
        T_diastolic = pulseTime(locs_peaks(2));
        xline(T_diastolic, 'k--', sprintf("%.2f s", T_diastolic), ...
            'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom', 'Color', [0.25 0.25 0.25], 'FontSize', 10);

        yline(peaks(2), 'k--', sprintf("%.1f %s", peaks(2), unit), ...
            'LineWidth', 1.5, 'LabelVerticalAlignment', 'top', 'Color', [0.25 0.25 0.25], 'FontSize', 10);
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
axis([axT(1), axT(2), axP(3) - 0.2 * (axP(4) - axP(3)), axP(4) + 0.1 * (axP(4) - axP(3))]);

ylabel(y_label);
xlabel('Time (s)');
pbaspect([1.618 1 1]);
set(gca, 'LineWidth', 2, 'Box', 'on');

% Save Results
exportgraphics(hFig, fullfile(ToolBox.path_png, ...
    sprintf("%s_ArterialWaveformAnalysis_%s.png", ToolBox.folder_name, name)), ...
    'Resolution', 300);

%

figure,
hold on
% plot(pulseTime, composite_signal, 'LineWidth', 2);

combined_composite = sum(composite_signal, 2);
combined_composite = rescale(combined_composite, min(signal, [], 'all'), max(signal, [], 'all')); % rescale at initial scale

plot(pulseTime, combined_composite, '--', 'LineWidth', 2);
plot(pulseTime, one_cycle_signal, 'k', 'LineWidth', 2);

% Configure axes
axis([0, pulseTime(end), axP(3) - 0.2 * (axP(4) - axP(3)), axP(4) + 0.1 * (axP(4) - axP(3))]);

ylabel(y_label);
xlabel('Time (s)');
pbaspect([1.618 1 1]);
set(gca, 'LineWidth', 2, 'Box', 'on');

% Save Results with harmonics
exportgraphics(hFig, fullfile(ToolBox.path_png, ...
    sprintf("%s_ArterialWaveformAnalysisHarmonics_%s.png", ToolBox.folder_name, name)), ...
    'Resolution', 300);

% Export to JSON (only for velocity signals)
if ~isBVR
    ToolBox.Output.add('SystoleDuration', dicroticNotchTime, 's');
    ToolBox.Output.add('DiastoleDuration', diastoleDuration, 's');
    ToolBox.Output.add('SystolicUpstroke', systolicUpstroke, unit);
    ToolBox.Output.add('SystolicDownstroke', systolicDownstroke, unit);
    ToolBox.Output.add('DiastolicRunoff', diastolicRunoff, unit);
    %     ToolBox.Output.add('DicroticNotchIndex', dicroticNotchIndex, unit);
end

% Close the figure if not needed
if ~strcmpi(get(hFig, 'Visible'), 'on')
    close(hFig);
end

end
