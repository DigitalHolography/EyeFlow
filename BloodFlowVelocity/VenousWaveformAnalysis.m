function VenousWaveformAnalysis(signal, t, sysIdxList, numInterp, name, ToolBox, cshiftn)
%VenousWaveformAnalysis Output figures and signal from venous analysis

if strcmp(name, "bvr")
    y_label = 'Blood Volume Rate (µL/min)';
    folder = 'crossSectionsAnalysis';
    unit = 'µL/min';
else
    y_label = 'Velocity (mm/s)';
    folder = 'bloodFlowVelocity';
    unit = 'mm/s';
end

% Cycle Analysis

try
    signal = double(wdenoise(signal));
catch
    signal = double(signal);
end

[one_cycle_signal, avgLength] = interpSignal(signal, sysIdxList, numInterp);

% we can use the cshiftn calculated for the veins to have a simultaneous interpolated waveforms plot
if nargin < 8
    [~, amin] = min(one_cycle_signal);
    cshiftn = mod(numInterp - amin, numInterp) + 1;
end

signal_shifted = circshift(one_cycle_signal, cshiftn);
signal_shifted(end + 1) = signal_shifted(1);

dt = (t(2) - t(1));
pulseTime = linspace(0, dt * avgLength, numInterp);
pulseTime(end + 1) = pulseTime(end) + dt;
figure("Visible", "off"), plot(pulseTime, signal_shifted, 'k', 'LineWidth', 2)
hold on

axis padded
axP = axis;
axis tight
axT = axis;
axis([axT(1), axT(2), 0, axP(4)])

ylabel(y_label)
xlabel('Time (s)')
pbaspect([1.618 1 1])
set(gca, 'LineWidth', 2), box on

[peaks, locs_peaks] = findpeaks(signal_shifted, 'MinPeakWidth', numInterp / 2);

Max = max(signal_shifted);
Min = min(signal_shifted);
Range = Max - Min;

% --- Ascent Time Detection (Fixed) ---
ascent_value = (Max - 0.05 * Range);

% Find zero-crossings where signal_shifted rises above ascent_value
crossings = find(diff(sign(signal_shifted - ascent_value)) > 0);

if ~isempty(crossings)
    index_ascent = crossings(1);
else
    % Fallback: Use the minimum if no crossing exists
    [~, index_ascent] = min(signal_shifted);
end

yline(ascent_value, 'k--', 'LineWidth', 2, 'LabelVerticalAlignment', 'bottom', 'Color', [0.4 0.4 0.4])
xline(pulseTime(index_ascent), 'k--', 'LineWidth', 2, 'LabelVerticalAlignment', 'bottom', 'Color', [0.4 0.4 0.4])

% --- Descent Time Detection (Fixed) ---
descent_value = (Min + 0.05 * Range);
% Find last crossing where signal_shifted falls below descent_value
crossings = find(diff(sign(signal_shifted - descent_value)) < 0);

if ~isempty(crossings)
    index_descent = crossings(end);
else
    % Fallback: Use the maximum if no crossing exists
    [~, index_descent] = max(signal_shifted);
end

yline(descent_value, 'k--', 'LineWidth', 2, 'LabelVerticalAlignment', 'bottom', 'Color', [0.4 0.4 0.4])
xline(pulseTime(index_descent), 'k--', 'LineWidth', 2, 'LabelVerticalAlignment', 'bottom', 'Color', [0.4 0.4 0.4])

scatter(pulseTime(locs_peaks), peaks, 'r')

if ~isempty(locs_peaks)
    T_peak = pulseTime(locs_peaks(1));
else
    T_peak = NaN;
end

T_ascent = pulseTime(index_ascent);
T_descent = pulseTime(index_descent);

exportgraphics(gca, fullfile(ToolBox.path_png, folder, sprintf("%s_VenousWaveformAnalysis_%s.png", ToolBox.folder_name, name)))

% Export to JSON
if ~strcmp(name, "bvr") % only for the velocity signal

    ToolBox.Outputs.add('TimetoPeakFromMinVein', T_peak, 's');
    ToolBox.Outputs.add('TimetoAscentFromMinVein', T_ascent, 's');
    ToolBox.Outputs.add('TimetoDescentToMinVein', T_descent, 's');

end

if strcmp(name, "bvr")
    % ToolBox.Outputs.add('SystoleDurationBvr', systoleDuration, 's'); % pour l'instant n'existe pas comme sortie car pas d'info en plus forcément
    % ToolBox.Outputs.add('DiastoleDurationBvr', diastoleDuration, 's');
    % ToolBox.Outputs.add('SystolicUpstrokeBvr', systolicUpstroke, unit);
    % ToolBox.Outputs.add('SystolicDownstrokeBvr', systolicDownstroke, unit);
    % ToolBox.Outputs.add('DiastolicRunoffBvr', diastolicRunoff, unit);
end

end
