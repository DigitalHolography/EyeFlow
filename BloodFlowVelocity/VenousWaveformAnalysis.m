function VenousWaveformAnalysis(signal, systolesIndexes, numInterp, name)
% VenousWaveformAnalysis Analyzes venous waveform signals (velocity or Flow Rate)
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

cDark = [0 0 1];
cLight = [0.5 0.5 1];

% Set parameters based on signal type
if strcmp(name, "bvr")
    y_label = 'Flow Rate (µL/min)';
    unit = 'µL/min';
    isBVR = true;
else
    y_label = 'Velocity (mm/s)';
    unit = 'mm/s';
    isBVR = false;
end

% Signal Preprocessing
try
    % Apply wavelet denoising if possible
    signal = double(wdenoise(signal, 4, ...
        'Wavelet', 'sym4', ...
        'DenoisingMethod', 'Bayes', ...
        'ThresholdRule', 'Median'));
catch
    signal = double(signal);
end

% Apply bandpass filter (0.5-15 Hz) as suggested
[b, a] = butter(4, 15 / (fs / 2), 'low');
filtered_signal = filtfilt(b, a, signal);

% Cycle Analysis
[one_cycle_signal, avgLength] = interpSignal(filtered_signal, systolesIndexes, numInterp);

% Create time vector for one cycle
dt = (t(2) - t(1));
pulseTime = linspace(0, dt * avgLength, numInterp);

% Feature Detection

% Find extremes
[peak, loc_peak] = max(one_cycle_signal);
[trough, loc_trough] = min(one_cycle_signal);

% Visualization
hFig = figure('Visible', 'off', 'Color', 'w');
plot(pulseTime, one_cycle_signal, 'k', 'LineWidth', 2)
hold on

% Add reference lines and annotations
T_peak = pulseTime(loc_peak);
xline(T_peak, 'k--', sprintf("%.2f s", T_peak), ...
    'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom', 'Color', [0.25 0.25 0.25], 'FontSize', 10);

xline(dt * avgLength, 'k--', sprintf("%.2f s", dt * avgLength), ...
    'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom', 'Color', [0.25 0.25 0.25], 'FontSize', 10);

yline(peak, 'k--', sprintf("%.1f %s", peak, unit), ...
    'LineWidth', 1.5, 'Color', [0.25 0.25 0.25], 'FontSize', 10);

T_Min = pulseTime(loc_trough);
xline(T_Min, 'k--', sprintf("%.2f s", T_Min), ...
    'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom', 'Color', [0.25 0.25 0.25], 'FontSize', 10);

yline(trough, 'b--', sprintf("%.1f %s", trough, unit), ...
    'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom', 'Color', [0.25 0.25 0.25], 'FontSize', 10);

padded_signal = [one_cycle_signal(floor(numInterp / 2) + 1:end), ...
                     one_cycle_signal, ...
                     one_cycle_signal(1:floor(numInterp / 2))];
padded_time = [linspace(- dt * avgLength / 2, -dt, round(numInterp / 2)), ...
                   pulseTime, ...
                   linspace(dt * (avgLength + 1), dt * (avgLength * 3/2), round(numInterp / 2))];

% Main signal and gradient plots
plot(padded_time, padded_signal, 'Color', [0.85 0.85 0.85], 'LineWidth', 2)
plot(pulseTime, one_cycle_signal, 'k', 'LineWidth', 2);

% Plot detected features
scatter(pulseTime(loc_peak), peak, 100, 'filled', 'MarkerFaceColor', cDark, 'MarkerEdgeColor', 'k');
scatter(pulseTime(loc_trough), trough, 100, 'filled', 'MarkerFaceColor', cLight, 'MarkerEdgeColor', 'k');

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
    sprintf("%s_VenousWaveformAnalysis_%s.png", ToolBox.folder_name, name)), ...
    'Resolution', 300);

% Export to JSON
if ~isBVR

    ToolBox.Outputs.add('TimetoPeakFromMinVein', T_peak-T_Min, 's');

end

if isBVR
    % ToolBox.Outputs.add('SystoleDurationBvr', systoleDuration, 's'); % pour l'instant n'existe pas comme sortie car pas d'info en plus forcément
    % ToolBox.Outputs.add('DiastoleDurationBvr', diastoleDuration, 's');
    % ToolBox.Outputs.add('SystolicUpstrokeBvr', systolicUpstroke, unit);
    % ToolBox.Outputs.add('SystolicDownstrokeBvr', systolicDownstroke, unit);
    % ToolBox.Outputs.add('DiastolicRunoffBvr', diastolicRunoff, unit);
end

% Close the figure if not needed
if ~strcmpi(get(hFig, 'Visible'), 'on')
    close(hFig);
end

end
