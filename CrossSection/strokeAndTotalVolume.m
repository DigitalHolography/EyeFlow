function strokeAndTotalVolume(mean_BvrT, mean_std_BvrT, systolesIndexes, numInterp, name)

ToolBox = getGlobalToolBox;
params = ToolBox.getParams;
saveFigures = params.saveFigures;
fs = 1 / (ToolBox.stride / ToolBox.fs / 1000);

% Apply low-pass filter (15 Hz)
[b, a] = butter(4, 15 / (fs / 2), 'low');
filtered_BvrT = filtfilt(b, a, mean_BvrT);

[interp_BvrT, avgLength, interp_std_BvrT] = interpSignal(filtered_BvrT, systolesIndexes, numInterp, mean_std_BvrT);

if isempty(interp_BvrT)
    return
end

dt = (avgLength / 60) / (numInterp - 1);
pulseTime = linspace(0, avgLength / 60, numInterp);
total_volume_value = sum(interp_BvrT) * dt / 60 * 1000; % in nL

[~, amin] = min(interp_BvrT);
[~, amax] = max(interp_BvrT);

if strcmp(name, 'artery')
    stroke_volume_value = sum(interp_BvrT(1:min(amax, numInterp))) * dt / 60 * 1000; % in nL
    stroke_volume_value = stroke_volume_value + (sum(interp_BvrT(max(amin, 1):end)) * dt / 60 * 1000); % in nL
else
    stroke_volume_value = sum(interp_BvrT(1:max(amin, 1))) * dt / 60 * 1000;
end

if contains(name, 'artery')
    ToolBox.Output.add('ArteryCycleVolume', total_volume_value, 'nL', h5path = '/Artery/FlowRate/CycleVolume/Total');
    ToolBox.Output.add('ArterySystolicFraction', stroke_volume_value / total_volume_value, h5path = '/Artery/FlowRate/CycleVolume/SystolicFraction');
    ToolBox.Output.add('ArteryDiastolicFraction', (1 - stroke_volume_value / total_volume_value), h5path = '/Artery/FlowRate/CycleVolume/DiastolicFraction');
elseif contains(name, 'vein')
    ToolBox.Output.add('VeinCycleVolume', total_volume_value, 'nL', h5path = '/Vein/FlowRate/CycleVolume/Total');
    ToolBox.Output.add('VeinSystolicFraction', stroke_volume_value / total_volume_value, h5path = '/Vein/FlowRate/CycleVolume/SystolicFraction');
    ToolBox.Output.add('VeinDiastolicFraction', (1 - stroke_volume_value / total_volume_value), h5path = '/Vein/FlowRate/CycleVolume/DiastolicFraction');
end

if saveFigures
    % Plotting
    f = figure("Visible", "off", "Color", "w");
    hold off

    % Retinal Stroke Volume
    hold on
    curve1 = interp_BvrT;
    curve2 = 0 * ones(size(curve1));
    ft2 = [pulseTime, fliplr(pulseTime)];
    inBetween = [curve1, fliplr(curve2)]';

    if strcmp(name, 'artery')
        cLight = [1, 1/2, 1/2];
        cDark = [1, 0, 0];
    else
        cLight = [1/2, 1/2, 1];
        cDark = [0, 0, 1];
    end

    fill(ft2, inBetween, cLight, 'EdgeColor', 'none');
    xline(pulseTime(end), '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 2)

    if strcmp(name, 'artery')
        % Remaining Stroke Volume
        hold on
        curve1 = interp_BvrT;
        curve1 = curve1(1:min(amax, numInterp));
        curve2 = 0 * ones(size(curve1));
        ft2 = [pulseTime(1:min(amax, numInterp)), fliplr(pulseTime(1:min(amax, numInterp)))];
        inBetween = [curve1, fliplr(curve2)]';
        xline(pulseTime(min(amax, numInterp)), '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 2)
        fill(ft2, inBetween, cDark, 'EdgeColor', 'none');

        % Remaining Stroke Volume
        curve1 = interp_BvrT;
        curve1 = curve1(max(amin, 1):end);
        curve2 = 0 * ones(size(curve1));
        ft2 = [pulseTime(max(amin, 1):end), fliplr(pulseTime(max(amin, 1):end))];
        inBetween = [curve1, fliplr(curve2)]';
        xline(pulseTime(max(amin, 1)), '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 2)
        fill(ft2, inBetween, cDark, 'EdgeColor', 'none');
    else
        % Remaining Stroke Volume
        hold on
        curve1 = interp_BvrT;
        curve1 = curve1(1:max(amin, 1));
        curve2 = 0 * ones(size(curve1));
        ft2 = [pulseTime(1:max(amin, 1)), fliplr(pulseTime(1:max(amin, 1)))];
        inBetween = [curve1, fliplr(curve2)]';
        xline(pulseTime(max(amin, 1)), '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 2)
        fill(ft2, inBetween, cDark, 'EdgeColor', 'none');
    end

    % Grey STD and Signal
    interp_BvrT2 = repmat(interp_BvrT, 1, 3);
    interp_std_BvrT2 = repmat(interp_std_BvrT, 1, 3);
    pulseTime2 = [flip(- pulseTime - dt), pulseTime, pulseTime(end) + dt + pulseTime];

    hold on
    curve1 = interp_BvrT2 + interp_std_BvrT2;
    curve2 = interp_BvrT2 - interp_std_BvrT2;
    ft2 = [pulseTime2, fliplr(pulseTime2)];
    inBetween = [curve1, fliplr(curve2)]';
    cSTD = [0.7 0.7 0.7];
    fill(ft2, inBetween, cSTD, 'EdgeColor', 'none');
    plot(pulseTime2, interp_BvrT2, '-k', 'LineWidth', 2);

    yline(0, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 2)
    xline(0, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 2)

    axis padded
    axP = axis;
    axis tight
    axT = axis;
    axis([axT(1), axT(2), -2, axP(4) * 1.2])
    lower_bound = pulseTime(1) -1/2 * pulseTime(end);
    upper_bound = 3/2 * pulseTime(end);
    xlim([lower_bound, upper_bound])

    box on
    set(gca, 'LineWidth', 2);
    set(gca, 'PlotBoxAspectRatio', [1.618, 1, 1])
    ax = gca;
    ax.LineStyleOrderIndex = 1; % Reset if needed
    ax.SortMethod = 'depth'; % Try changing sorting method
    ax.Layer = 'top'; % This may help in some cases

    ylabel('Flow Rate (µL/min)')
    xlabel('Time (s)')

    dim = [0.2 0.5 0.3 0.3];

    str = sprintf("Retinal Stroke Volume : %02.0f nL\nTotal Volume : %02.0f nL", ...
        stroke_volume_value, total_volume_value);
    annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on', ...
        'EdgeColor', 'none', 'BackgroundColor', 'w');

    exportgraphics(gca, fullfile(ToolBox.path_png, sprintf("%s_strokeAndTotalVolume_%s.png", ToolBox.folder_name, name)))
    exportgraphics(gca, fullfile(ToolBox.path_eps, sprintf("%s_strokeAndTotalVolume_%s.eps", ToolBox.folder_name, name)))

    close(f);
end

end
