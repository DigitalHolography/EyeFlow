function ArterialWaveformAnalysis(signal, signal_ste, t, systolesIndexes, numInterp, name, ToolBox)

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

signal = double(signal);
[one_cycle_signal, avgLength, ~] = interpSignal(signal, systolesIndexes, numInterp, signal_ste);

[~, amin] = min(one_cycle_signal);
cshiftn = mod(numInterp - amin, numInterp) + 1;
signal_shifted = circshift(one_cycle_signal, cshiftn);
signal_shifted(end + 1) = signal_shifted(1);

dt = (t(2) - t(1));
pulseTime = linspace(0, dt * avgLength, numInterp);
pulseTime(end + 1) = pulseTime(end) + dt;

figure, plot(pulseTime, signal_shifted, 'k', 'LineWidth', 2)
hold on

axis padded
axP = axis;
axis tight
axT = axis;
axis([axT(1), axT(2), 0, axP(4)])

ylabel(y_label)
xlabel('Time (s)')
pbaspect([1.618 1 1])
fontsize(gca, 14, 'points')
set(gca, 'LineWidth', 2), box on

min_peak_height = max(signal_shifted) * 0.3; % Adaptive threshold
min_peak_distance = floor(length(signal_shifted) / 4); % Minimum distance between peaks

f = fit(pulseTime.', signal_shifted.', 'gauss4');
y = feval(f, pulseTime);
error = rmse(rescale(signal_shifted), rescale(y'));
[peaks, locs_peaks] = findpeaks(y, 'MinPeakHeight', min_peak_height, 'MinPeakDistance', min_peak_distance);
[notch, locs_notch] = min(signal_shifted(locs_peaks(1):locs_peaks(2)));

scatter(pulseTime(locs_peaks), peaks, 'r')

if length(peaks) > 1
    [notch, locs_notch] = min(signal_shifted(locs_peaks(1):locs_peaks(2)));
    locs_notch = locs_notch + locs_peaks(1);
    scatter(pulseTime(locs_notch), notch, 'b')
end

% T_pulse = dt * avgLength;
T_peak = pulseTime(locs_peaks(1));

xline(T_peak, 'k--', sprintf("Tpeak = %.2f s", T_peak), 'LineWidth', 2, 'LabelVerticalAlignment', 'bottom', 'Color', [0.4 0.4 0.4]);
yline(peaks(1), 'k--', sprintf("Peak Systolic = %.2f %s", peaks(1), unit), 'LineWidth', 2, 'LabelVerticalAlignment', 'bottom', 'Color', [0.4 0.4 0.4]);
systolicUpstroke = peaks(1) - signal_shifted(1);

if length(locs_peaks) > 1
    T_notch_end = pulseTime(locs_peaks(2));
    xline(T_notch_end, 'k--', sprintf("Tnotch %.2f s", T_notch_end), 'LineWidth', 2, 'LabelVerticalAlignment', 'bottom', 'Color', [0.4 0.4 0.4]);
    yline(peaks(2), 'k--', sprintf("Peak Diastolic %.2f %s", peaks(2), unit), 'LineWidth', 2, 'LabelVerticalAlignment', 'bottom', 'Color', [0.4 0.4 0.4]);

    T_notch = pulseTime(locs_notch);
    xline(T_notch, 'k--', sprintf("Systolic Phase Duration %.2f s", T_notch), 'LineWidth', 2, 'LabelVerticalAlignment', 'bottom', 'Color', [0.4 0.4 0.4]);
    yline(notch(1), 'k--', sprintf("Dicrotic Notch %.2f %s", notch(1), unit), 'LineWidth', 2, 'LabelVerticalAlignment', 'bottom', 'Color', [0.4 0.4 0.4]);

    systoleDuration = T_notch - pulseTime(1);
    diastoleDuration = pulseTime(end - 1) - T_notch;
    systolicDownstroke = peaks(1) - notch(1);
    diastolicRunoff = notch(1) - signal_shifted(1);
end

exportgraphics(gca, fullfile(ToolBox.path_png, folder, sprintf("%s_ArterialWaveformAnalysis_%s.png", ToolBox.main_foldername, name)))

% Export to JSON

ToolBox.Outputs.add('SystoleDuration', systoleDuration, 's');
ToolBox.Outputs.add('DiastoleDuration', diastoleDuration, 's');
ToolBox.Outputs.add('SystolicUpstroke', systolicUpstroke, unit);
ToolBox.Outputs.add('SystolicDownstroke', systolicDownstroke, unit);
ToolBox.Outputs.add('DiastolicRunoff', diastolicRunoff, unit);

figure, plot(f, pulseTime, signal_shifted)

end
