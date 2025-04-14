function ArterialWaveformAnalysis(signal, signal_ste, t, systolesIndexes, numInterp, name, ToolBox)

if strcmp(name, "bvr")
    y_label = 'Blood Volume Rate (µL/min)';
else
    y_label = 'Velocity (mm/s)';
end

% Cycle Analysis

[one_cycle_signal, avgLength, ~] = interpSignal(signal, systolesIndexes, numInterp, signal_ste);

[~, amin] = min(one_cycle_signal);
cshiftn = mod(numInterp - amin, numInterp) + 1;
signal_shifted = circshift(one_cycle_signal, cshiftn);
signal_shifted(end + 1) = signal_shifted(1);

dt = (t(2) - t(1));
pulseTime = linspace(0, dt * avgLength, numInterp);
pulseTime(end + 1) = pulseTime(end) + dt;

figure("Visible", 'off'), plot(pulseTime, signal_shifted, 'k', 'LineWidth', 2)
hold on

axis padded
axP = axis;
axis tight
axT = axis;
axis([axT(1), axT(2), 0, axP(4)])

ylabel('Blood Volume Rate (µL/min)')
xlabel('Time (s)')
pbaspect([1.618 1 1])
fontsize(gca, 14, 'points')
set(gca, 'LineWidth', 2), box on

min_peak_height = max(signal_shifted) * 0.3; % Adaptive threshold
min_peak_distance = floor(length(signal_shifted) / 3); % Minimum distance between peaks

[peaks, locs_peaks] = findpeaks(signal_shifted, 'MinPeakHeight', min_peak_height, 'MinPeakDistance', min_peak_distance);
scatter(pulseTime(locs_peaks), peaks, 'r')
if length(peaks) > 1
    [notch, locs_notch] = min(signal_shifted(locs_peaks(1):locs_peaks(2)));
    locs_notch = locs_notch + locs_peaks(1);
    scatter(pulseTime(locs_notch), notch, 'b')
end

T_pulse = dt * avgLength;
T_peak = pulseTime(locs_peaks(1));

xline(T_peak, 'k--', sprintf("Tpeak = %.2f s", T_peak), 'LineWidth', 2, 'LabelVerticalAlignment', 'bottom');
yline(peaks(1), 'k--', sprintf("Peak Systolic = %.2f µL/min", peaks(1)), 'LineWidth', 2, 'LabelVerticalAlignment', 'bottom');
systolicUpstroke = peaks(1) - signal_shifted(1);

if length(locs_peaks) > 1
    T_notch_end = pulseTime(locs_peaks(2));
    xline(T_notch_end, 'k--', sprintf("Tnotch %.2f s", T_notch_end), 'LineWidth', 2, 'LabelVerticalAlignment', 'bottom');
    yline(peaks(2), 'k--', sprintf("Peak Diastolic %.2f µL/min", peaks(2)), 'LineWidth', 2, 'LabelVerticalAlignment', 'bottom');

    T_notch = pulseTime(locs_notch);
    xline(T_notch, 'k--', sprintf("Systolic Phase Duration %.2f s", T_notch), 'LineWidth', 2, 'LabelVerticalAlignment', 'bottom');
    yline(notch(1), 'k--', sprintf("Dicrotic Notch %.2f µL/min", notch(1)), 'LineWidth', 2, 'LabelVerticalAlignment', 'bottom');

    systoleDuration = T_notch - pulseTime(1);
    diastoleDuration = pulseTime(end - 1) - T_notch;
    systolicDownstroke = peaks(1) - notch(1);
    diastolicRunoff = notch(1) - signal_shifted(1);
end

exportgraphics(gca, fullfile(ToolBox.path_png, 'crossSectionsAnalysis', sprintf("%s_ArterialWaveformAnalysis_%s.png", ToolBox.main_foldername, name)))

end
