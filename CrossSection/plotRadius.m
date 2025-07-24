function [Q_t, dQ_t] = plotRadius(radiusQ, radiusQSE, fullTime, idx_start, idx_end, name)
% Get global toolbox and parameters
ToolBox = getGlobalToolBox;
params = ToolBox.getParams;
numCircles = params.json.CrossSectionsAnalysis.NumberOfCircles; % Number of circles (radii)
r1 = params.json.SizeOfField.SmallRadiusRatio;
r2 = params.json.SizeOfField.BigRadiusRatio;
dr = (r2 - r1) / numCircles;
rad = linspace(r1, r2 - dr, numCircles);

% Define color for shaded regions
Color_std = [0.7 0.7 0.7];

% Compute time-averaged mean and standard deviation
Q_r = squeeze(mean(radiusQ(:, idx_start:idx_end), 2))'; % Mean over time
N = length(idx_start:idx_end);
dQ_r = squeeze(sqrt(sum(radiusQSE(:, idx_start:idx_end) .^ 2, 2)))' / N; % RMS of uncertainties

% Create shaded region for uncertainty
curve1 = Q_r + dQ_r; % Upper bound
curve2 = Q_r - dQ_r; % Lower bound
rad2 = [rad, fliplr(rad)]; % X-values for fill
inBetween = [curve1, fliplr(curve2)]; % Y-values for fill

% Plot time-averaged Volume Rate vs. radius
figure("Visible", "off");
fill(rad2, inBetween, Color_std, 'EdgeColor', 'none'); % Shaded region
hold on;
plot(rad, curve1, "Color", Color_std, 'LineWidth', 2); % Upper bound
plot(rad, curve2, "Color", Color_std, 'LineWidth', 2); % Lower bound
plot(rad, Q_r, '-k', 'LineWidth', 2); % Mean curve
yline(mean(Q_r), '--k', 'LineWidth', 2, 'label', sprintf("%0.2f µL/min", mean(Q_r))); % Horizontal mean line

% Format plot
axis padded;
axP = axis;
axis tight;
axT = axis;
axis([axT(1), axT(2), 0, 1.07 * axP(4)]);

box on
set(gca, 'LineWidth', 2);
set(gca, 'PlotBoxAspectRatio', [2.5, 1, 1])
ax = gca;
ax.LineStyleOrderIndex = 1; % Reset if needed
ax.SortMethod = 'depth'; % Try changing sorting method
ax.Layer = 'top'; % This may help in some cases

ylabel('Volume Rate (µL/min)');
xlabel('Radius (pixels)');
title("Time-Averaged Volume Rate");

% Export plot
exportgraphics(gca, fullfile(ToolBox.path_png, sprintf("%s_mean_%s_radius.png", ToolBox.folder_name, name)));

% Plot radial variations of Volume Rate over time
figure("Visible", "off");
hold on;

for circleIdx = 1:numCircles
    plot(fullTime, radiusQ(circleIdx, :));
end

plot(fullTime, mean(radiusQ, 1), 'k', 'LineWidth', 2);

axis padded;
axP = axis;
axis tight;
axT = axis;
axis([axT(1), axT(2), 0, 1.07 * axP(4)]);
box on;
ylabel('Volume Rate (µL/min)');
xlabel('Time (s)');
title("Radial Variations of Volume Rate");
set(gca, 'PlotBoxAspectRatio', [2.5 1 1], 'LineWidth', 2);

% Export plot
exportgraphics(gca, fullfile(ToolBox.path_png, sprintf("%s_variance_%s_time.png", ToolBox.folder_name, name)));

% Compute total Volume Rate over circles
Q_t = squeeze(mean(radiusQ, 1)); % Mean over circles
N = size(radiusQSE, 1);
dQ_t = squeeze(sqrt(sum(radiusQSE .^ 2, 1))) / N; % RMS of uncertainties

% Compute statistics for the time range
mean_Q = mean(Q_t(idx_start:idx_end)); % Time-averaged mean
[max_Q, amax] = max(Q_t(idx_start:idx_end)); % Maximum value in the range
N = length(idx_start:idx_end);
mean_dQ = sqrt(sum(dQ_t(idx_start:idx_end) .^ 2)) / N; % RMS of the uncertainty
[min_Q, amin] = min(Q_t(idx_start:idx_end)); % Minimum value in the range
max_dQ = dQ_t(amax);
min_dQ = dQ_t(amin);

% Plot total Volume Rate over time
figure("Visible", "off");
curve1 = Q_t + dQ_t; % Upper bound
curve2 = Q_t - dQ_t; % Lower bound
ft2 = [fullTime, fliplr(fullTime)]; % X-values for fill
inBetween = [curve1, fliplr(curve2)]; % Y-values for fill

hold on;
fill(ft2, inBetween, Color_std, 'EdgeColor', 'none'); % Shaded region
yline(0, 'k-', 'LineWidth', 2); % Zero line
plot(fullTime, curve1, "Color", Color_std, 'LineWidth', 2); % Upper bound
plot(fullTime, curve2, "Color", Color_std, 'LineWidth', 2); % Lower bound
plot(fullTime, Q_t, '-k', 'LineWidth', 2); % Mean curve
yline(mean_Q, '--k', 'LineWidth', 2); % Horizontal mean line

% Mark the time range used for averaging
plot(fullTime(idx_start), 1.07 * max_Q, 'k|', 'MarkerSize', 10, 'LineWidth', 2);
plot(fullTime(idx_end), 1.07 * max_Q, 'k|', 'MarkerSize', 10, 'LineWidth', 2);
plot(fullTime(idx_start:idx_end), repmat(1.07 * max_Q, 1, idx_end - idx_start + 1), '-k', 'LineWidth', 2);

% Format plot
axis padded;
axP = axis;
axis tight;
axT = axis;
axis([axT(1), axT(2), axP(3), 1.07 * axP(4)]);

box on
set(gca, 'LineWidth', 2);
set(gca, 'PlotBoxAspectRatio', [2.5, 1, 1])
ax = gca;
ax.LineStyleOrderIndex = 1; % Reset if needed
ax.SortMethod = 'depth'; % Try changing sorting method
ax.Layer = 'top'; % This may help in some cases

ylabel('Volume Rate (µL/min)');
xlabel('Time (s)');
title(sprintf("Total Volume Rate (Avg. %0.2f µL/min)", mean_Q));

% Export plot
exportgraphics(gca, fullfile(ToolBox.path_png, sprintf("%s_allrad_%s_time.png", ToolBox.folder_name, name)));
exportgraphics(gca, fullfile(ToolBox.path_eps, sprintf("%s_allrad_%s_time.eps", ToolBox.folder_name, name)));

% Write results to a text file
fileID = fopen(fullfile(ToolBox.path_txt, strcat(ToolBox.folder_name, '_', 'EF_main_outputs', '.txt')), 'a');
fprintf(fileID, 'Flow Rate %s : %f (µL/min) \r\n', name, mean_Q);
fprintf(fileID, 'Flow Rate Standard Deviation %s : %f (µL/min) \r\n', name, mean_dQ);
fclose(fileID);

% Write results to json
ToolBox.outputs.(sprintf('FlowRate%s', name)) = mean_Q;
ToolBox.outputs.(sprintf('FlowRateStd%s', name)) = mean_dQ;

% New
if contains(name, 'Vein')
    ToolBox.Outputs.add('VenousMeanVolumeRate', mean_Q, 'µL/min', mean_dQ);
    ToolBox.Outputs.add('VenousMaximumVolumeRate', max_Q, 'µL/min', max_dQ);
    ToolBox.Outputs.add('VenousMinimumVolumeRate', min_Q, 'µL/min', min_dQ);
elseif contains(name, 'Artery')
    ToolBox.Outputs.add('ArterialMeanVolumeRate', mean_Q, 'µL/min', mean_dQ);
    ToolBox.Outputs.add('ArterialMinimumVolumeRate', min_Q, 'µL/min', min_dQ);
    ToolBox.Outputs.add('ArterialMaximumVolumeRate', max_Q, 'µL/min', max_dQ);
end

end
