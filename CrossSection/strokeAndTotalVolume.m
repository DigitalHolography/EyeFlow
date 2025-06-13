function strokeAndTotalVolume(mean_BvrT, mean_std_BvrT, systolesIndexes, fullTime, numInterp, name)

ToolBox = getGlobalToolBox;

figAspect("Visible", 1);

[interp_BvrT, avgLength, interp_std_BvrT] = interpSignal(mean_BvrT, systolesIndexes, numInterp, mean_std_BvrT);

if isempty(interp_BvrT)
    return
end

dt = (avgLength / 60)/(numInterp - 1);
pulseTime = linspace(0, avgLength / 60, numInterp);

[mindiastole_bvr_value, amin] = min(interp_BvrT);
[maxsystole_bvr_value, amax] = max(interp_BvrT);

hold off

% Retinal Stroke Volume
hold on
curve1 = interp_BvrT;
curve2 = 0 * ones(size(curve1));
ft2 = [pulseTime, fliplr(pulseTime)];
inBetween = [curve1, fliplr(curve2)]';

if strcmp(name, 'Artery')
    cLight = [1, 1/2, 1/2];
    cDark = [1, 0, 0];
else
    cLight = [0, 0, 1];
    cDark = [1/2, 1/2, 1];
end

fill(ft2, inBetween, cLight, 'EdgeColor', 'none');
xline(pulseTime(end), '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 2)

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
xline(0, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 2)

axis padded
axP = axis;
axis tight
axT = axis;
axis([axT(1), axT(2), -2, axP(4) * 1.2])
lower_bound = pulseTime(1) -1/2 * pulseTime(end);
upper_bound = 3/2 * pulseTime(end);
xlim([lower_bound, upper_bound])

ylabel('Blood Volume Rate (µL/min)')
xlabel('Time (s)')
stroke_volume_value = sum(interp_BvrT(1:min(amax, numInterp))) * dt / 60 * 1000; % in nL
stroke_volume_value = stroke_volume_value + (sum(interp_BvrT(max(amin, 1):end)) * dt / 60 * 1000); % in nL
total_volume_value = sum(interp_BvrT) * dt / 60 * 1000;

dim = [0.2 0.5 0.3 0.3];
if strcmp(name, 'Artery')
    str = sprintf("Retinal Stroke Volume : %02.0f nL and Total Volume : %02.0f nL", stroke_volume_value, total_volume_value);
    annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on', 'BackgroundColor', 'w');
else
    str = sprintf("Total Volume : %02.0f nL", total_volume_value);
    annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on', 'BackgroundColor', 'w');
end

exportgraphics(gca, fullfile(ToolBox.path_png, 'crossSectionsAnalysis', ...
    sprintf("%s_strokeAndTotalVolume_%s.png", ToolBox.folder_name, name)))
exportgraphics(gca, fullfile(ToolBox.path_eps, 'crossSectionsAnalysis', ...
    sprintf("%s_strokeAndTotalVolume_%s.eps", ToolBox.folder_name, name)))

fileID = fopen(fullfile(ToolBox.path_txt, sprintf('%s_EF_main_outputs.txt', ToolBox.folder_name)), 'a');
fprintf(fileID, 'MaxSystole Blood Volume Rate Artery : %f (µL/min) \r\n', maxsystole_bvr_value);
fprintf(fileID, 'MinDiastole Blood Volume Rate Artery : %f (µL/min) \r\n', mindiastole_bvr_value);
fprintf(fileID, 'Stroke Volume Artery : %f (nL) \r\n', stroke_volume_value);
fprintf(fileID, 'Total Volume Artery : %f (nL) \r\n', total_volume_value);
fclose(fileID);

ToolBox.outputs.(sprintf('MaxSystoleFlowRateArtery')) = maxsystole_bvr_value;
ToolBox.outputs.(sprintf('MinDiastoleFlowRateArtery')) = mindiastole_bvr_value;
ToolBox.outputs.(sprintf('StrokeVolumeArteryArtery')) = stroke_volume_value;
ToolBox.outputs.(sprintf('TotalVolumeArteryArtery')) = total_volume_value;

if contains(name, 'Artery')
    ToolBox.Outputs.add('ArterialCycleVolume', total_volume_value, 'nL');
    ToolBox.Outputs.add('ArterialSystolicFraction', stroke_volume_value / total_volume_value, '');
    ToolBox.Outputs.add('ArterialDiastolicFraction', (1 - stroke_volume_value / total_volume_value), '');
elseif contains(name, 'Vein')
    ToolBox.Outputs.add('VenousCycleVolume', total_volume_value, 'nL');
    ToolBox.Outputs.add('VenousSystolicFraction', stroke_volume_value / total_volume_value, '');
    ToolBox.Outputs.add('VenousDiastolicFraction', (1 - stroke_volume_value / total_volume_value), '');
end

close all

end
