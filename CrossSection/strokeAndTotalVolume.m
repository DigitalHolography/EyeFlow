function strokeAndTotalVolume(mean_BvrT, mean_std_BvrT, systolesIndexes, fullTime, numInterp, name)

ToolBox = getGlobalToolBox;

figAspect;

[interp_BvrT, avgLength, interp_std_BvrT] = interpSignal(mean_BvrT, systolesIndexes, numInterp, mean_std_BvrT);

if isempty(interp_BvrT)
    return
end

dt = (fullTime(2) - fullTime(1));
pulseTime = dt * (1:numInterp) * avgLength / numInterp;

[mindiastole_bvr_value, amin] = min(interp_BvrT);
[maxsystole_bvr_value, amax] = max(interp_BvrT);
cshiftn = mod(numInterp - amin + 1, numInterp);

hold off

% Retinal Stroke Volume
hold on
curve1 = circshift(interp_BvrT, cshiftn);
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
xline(pulseTime(end), 'k--', 'LineWidth', 2)

% Remaining Stroke Volume
hold on
curve1 = circshift(interp_BvrT, cshiftn);
curve1 = curve1(1:min(amax + cshiftn, numInterp));
curve2 = 0 * ones(size(curve1));
ft2 = [pulseTime(1:min(amax + cshiftn, numInterp)), fliplr(pulseTime(1:min(amax + cshiftn, numInterp)))];
inBetween = [curve1, fliplr(curve2)]';
xline(pulseTime(min(amax + cshiftn, numInterp)), 'k--', 'LineWidth', 2)
fill(ft2, inBetween, cDark, 'EdgeColor', 'none');

% Grey STD and Signal
interp_BvrT2 = repmat(interp_BvrT, 1, 3);
interp_std_BvrT2 = repmat(interp_std_BvrT, 1, 3);
pulseTime2 = dt * (-numInterp + 1:numInterp * 2) * avgLength / numInterp;

hold on
curve1 = circshift(interp_BvrT2, cshiftn) + circshift(interp_std_BvrT2, cshiftn);
curve2 = circshift(interp_BvrT2, cshiftn) - circshift(interp_std_BvrT2, cshiftn);
ft2 = [pulseTime2, fliplr(pulseTime2)];
inBetween = [curve1, fliplr(curve2)]';
cSTD = [0.7 0.7 0.7];
fill(ft2, inBetween, cSTD, 'EdgeColor', 'none');
plot(pulseTime2, circshift(interp_BvrT2, cshiftn), '-k', 'LineWidth', 2);

yline(0, 'k--', 'LineWidth', 2)
xline(0, 'k--', 'LineWidth', 2)
xline(0, 'k--', 'LineWidth', 2)

axis padded
axP = axis;
axis tight
axT = axis;
axis([axT(1), axT(2), -2, axP(4) * 1.07])
lower_bound = pulseTime(1) -1/2 * pulseTime(end);
upper_bound = 3/2 * pulseTime(end);
xlim([lower_bound, upper_bound])

ylabel('Blood Volume Rate (µL/min)')
xlabel('Time (s)')
ccinterpBvrT = circshift(interp_BvrT, cshiftn);
dt2 = pulseTime2(2) - pulseTime2(1);
stroke_volume_value = sum(ccinterpBvrT(1:min(amax + cshiftn, numInterp))) * dt2 / 60 * 1000; % in nL
total_volume_value = sum(ccinterpBvrT) * dt2 / 60 * 1000;

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

% ArterialWaveformAnalysis(mean_BvrT, mean_std_BvrT, fullTime, systolesIndexes, numInterp, sprintf('%s_bvr', name), ToolBox)

end
