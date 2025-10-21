function [sys_idx_list, pulse_artery_filtered, sys_max_list, sys_min_list] = find_systole_index(pulse_artery, options)
% FIND_SYSTOLE_INDEX Identifies systole peaks in the pulse signal.
% Inputs:
%   - pulseArtery: The extracted pulse signal
%   - options: Struct containing optional parameters
%       - pulseVein: (Optional) The extracted vein pulse signal
%       - savepng: (Optional) Boolean flag to save diagnostic plots
% Outputs:
%   - sys_idx_list: Indices of systole peaks
%   - pulse_artery_filtered: The extracted pulse signal
%   - sys_max_list: Indices of local maxima within each cycle
%   - sys_min_list: Indices of local minima within each cycle

arguments
    pulse_artery (:, 1) double
    options.pulseVein (:, 1) double = []
    options.savepng (1, 1) logical = true
    options.lowpass_freq (1, 1) double = 15
end

% Input validation
ToolBox = getGlobalToolBox();
fs = ToolBox.fs * 1000 / ToolBox.stride; % Convert to Hz

flagVein = ~isempty(options.pulseVein);

% Step 1: Extract pulse signal
[b, a] = butter(4, options.lowpass_freq / (fs / 2), 'low');
pulse_artery_filtered = filtfilt(b, a, pulse_artery);

if flagVein
    pulse_vein_filtered = filtfilt(b, a, options.pulseVein);
end

% Step 2: Compute derivative
diff_artery_signal = gradient(pulse_artery_filtered);

if flagVein
    diff_vein_signal = gradient(pulse_vein_filtered);
end

% Step 3: Detect peaks
min_peak_height = prctile(diff_artery_signal, 95);
min_peak_distance = floor(length(pulse_artery_filtered) / 7); % Minimum distance between peaks (0.5 seconds for a 3.5s video)
[~, sys_idx_list] = findpeaks(diff_artery_signal, 'MinPeakHeight', min_peak_height, 'MinPeakDistance', min_peak_distance);

% Step 4: Validate peaks
sys_idx_list = validate_peaks(sys_idx_list, 10);

% Step 5: Find local maxima and minima
num_peaks = numel(sys_idx_list);
sys_max_list = zeros(num_peaks, 1);
sys_min_list = zeros(num_peaks, 1);

for i = 1:(numel(sys_idx_list) - 1)
    L = sys_idx_list(i + 1) - sys_idx_list(i);
    D = round(L / 2);
    % Find the maximum within the current cycle
    [~, amax] = max(pulse_artery_filtered(sys_idx_list(i):sys_idx_list(i) + D));
    sys_max_list(i) = sys_idx_list(i) + amax - 1;

    % Find the minimum within the current cycle
    [~, amin] = min(pulse_artery_filtered(sys_idx_list(i) + D:sys_idx_list(i + 1) - 1));
    sys_min_list(i + 1) = sys_idx_list(i) + amin - 1 + D;
end

% Find the minimum before the first cycle
[~, amin] = min(pulse_artery_filtered(1:sys_idx_list(1)));
sys_min_list(1) = amin;

% Find the maximum after the end cycle
[~, amax] = max(pulse_artery_filtered(sys_idx_list(i + 1):end));
sys_max_list(i + 1) = sys_idx_list(i + 1) + amax - 1;

sys_max_list = sys_max_list';
sys_min_list = sys_min_list';

% Step 6: Error handling
if isempty(sys_idx_list)
    error('No systole peaks detected. Check signal quality or adjust parameters.');
end

% % FOR DEBUG
if options.savepng

    ToolBox = getGlobalToolBox();
    T = ToolBox.stride / ToolBox.fs / 1000;
    t = ToolBox.Cache.t;

    figure(Visible = 'off');
    hold on
    plot(t, diff_artery_signal, 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5)
    plot(t, pulse_artery_filtered, 'k-', 'LineWidth', 1.5);
    scatter((sys_max_list - 1) * T, pulse_artery_filtered(sys_max_list), 'r', "filled")
    scatter((sys_min_list - 1) * T, pulse_artery_filtered(sys_min_list), 'b', "filled")
    scatter((sys_idx_list - 1) * T, pulse_artery_filtered(sys_idx_list), 'k', "filled")

    xline((sys_idx_list - 1) * T, 'k--', 'LineWidth', 1.5)
    xline((sys_min_list - 1) * T, 'b--')
    xline((sys_max_list - 1) * T, 'r--')
    hold off

    axis padded;
    axP = axis;
    axis tight;
    axT = axis;
    axis([axT(1), axT(2), axP(3), axP(4)]);
    box on
    set(gca, 'LineWidth', 2, 'PlotBoxAspect', [1.618 1 1])
    xlabel("Time (s)")
    ylabel("Velocity (mm/s)")

    exportgraphics(gca, fullfile(ToolBox.path_png, sprintf("%s_%s", ToolBox.folder_name, 'find_systoles_indices_artery.png')))

    if flagVein
        figure(Visible = 'off');
        hold on
        plot(t, diff_vein_signal, 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5)
        plot(t, pulse_vein_filtered, 'k-', 'LineWidth', 1.5);
        scatter((sys_max_list - 1) * T, pulse_vein_filtered(sys_max_list), 'r', "filled")
        scatter((sys_min_list - 1) * T, pulse_vein_filtered(sys_min_list), 'b', "filled")
        scatter((sys_idx_list - 1) * T, pulse_vein_filtered(sys_idx_list), 'k', "filled")

        xline((sys_idx_list - 1) * T, 'k--', 'LineWidth', 1.5)
        xline((sys_min_list - 1) * T, 'b--')
        xline((sys_max_list - 1) * T, 'r--')
        hold off

        axis padded;
        axP = axis;
        axis tight;
        axT = axis;
        axis([axT(1), axT(2), axP(3), axP(4)]);
        box on
        set(gca, 'LineWidth', 2, 'PlotBoxAspect', [1.618 1 1])
        xlabel("Time (s)")
        ylabel("Velocity (mm/s)")

        exportgraphics(gca, fullfile(ToolBox.path_png, sprintf("%s_%s", ToolBox.folder_name, 'find_systoles_indices_vein.png')))
    end

end

end

%% **Validate Peaks (Removes peaks that are too close)**
function sys_idx_list = validate_peaks(sys_idx_list, min_distance)
i = 1;

while i < numel(sys_idx_list)

    if sys_idx_list(i + 1) - sys_idx_list(i) < min_distance
        sys_idx_list(i + 1) = [];
    else
        i = i + 1;
    end

end

end
