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
t = linspace(0, numFrames * ToolBox.stride / ToolBox.fs / 1000, numFrames);

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
hFig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1000 620]);
hold on;

% Add reference lines and annotations
T_peak = pulseTime(locs_peaks(1));
xline(T_peak, 'k--', sprintf("Max = %.2f s", T_peak), ...
    'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom', 'Color', [0.25 0.25 0.25], 'FontSize', 10);

xline(dt * avgLength, 'k--', sprintf("Total = %.2f s", dt * avgLength), ...
    'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom', 'Color', [0.25 0.25 0.25], 'FontSize', 10);

yline(peaks(1), 'k--', sprintf("Peak Systolic = %.1f %s", peaks(1), unit), ...
    'LineWidth', 1.5, 'Color', [0.25 0.25 0.25], 'FontSize', 10);

T_Min = pulseTime(endMinLoc);
xline(T_Min, 'k--', sprintf("Min = %.2f s", T_Min), ...
    'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom', 'Color', [0.25 0.25 0.25], 'FontSize', 10);

yline(endMinVal, 'b--', sprintf("End Minimum = %.1f %s", endMinVal, unit), ...
    'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom', 'Color', [0.25 0.25 0.25], 'FontSize', 10);

if ~isnan(notch)
    T_notch = pulseTime(locs_notch);
    xline(T_notch, 'k--', sprintf("Notch = %.2f s", T_notch), ...
        'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom', 'Color', [0.25 0.25 0.25], 'FontSize', 10);

    yline(notch, 'k--', sprintf("Notch = %.1f %s", notch, unit), ...
        'LineWidth', 1.5, 'Color', [0.25 0.25 0.25], 'FontSize', 10);

    if length(locs_peaks) > 1
        T_diastolic = pulseTime(locs_peaks(2));
        xline(T_diastolic, 'k--', sprintf("Max2 = %.2f s", T_diastolic), ...
            'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom', 'Color', [0.25 0.25 0.25], 'FontSize', 10);

        yline(peaks(2), 'k--', sprintf("Peak Diastolic = %.1f %s", peaks(2), unit), ...
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
scatter(pulseTime(locs_peaks), peaks, 100, 'r', 'filled', 'MarkerEdgeColor', 'k');
scatter(pulseTime(endMinLoc), endMinVal, 100, 'b', 'filled', 'MarkerEdgeColor', 'k');

if ~isnan(notch)
    scatter(pulseTime(locs_notch), notch, 100, 'b', 'filled', 'MarkerEdgeColor', 'k');
end

% Configure axes
axis tight;
axT = axis;
axis([axT(1), axT(2), axT(3) - 0.1 * (axT(4) - axT(3)), axT(4) + 0.1 * (axT(4) - axT(3))]);

ylabel(y_label);
xlabel('Time (s)');
pbaspect([1.618 1 1]);
set(gca, 'LineWidth', 2, 'FontSize', 24, 'Box', 'on');

%% Save Results
try
    exportgraphics(hFig, fullfile(ToolBox.path_png, folder, ...
        sprintf("%s_ArterialWaveformAnalysis_%s.png", ToolBox.folder_name, name)), ...
        'Resolution', 300);
catch ME
    warning('FailedToSaveFigure', 'Could not save figure: %s', ME.message);
end

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
