function [U_t, U_t_SE] = plotRadius(radius_U, radius_U_SE, fullTime, idx_start, idx_end, vessel_name, type_name)
% plotRadius - Plots the average volume rate per radius and its uncertainty
%   This function computes and plots the average volume rate and its
%   uncertainty over a specified time range for a given radius.

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

if strcmp(type_name, 'velocity')
    y_label = 'Velocity (mm/s)';
else
    y_label = 'Volume Rate (µL/min)';
end

% Compute time-averaged mean and standard deviation
U_r = squeeze(mean(radius_U(:, idx_start:idx_end), 2))'; % Mean over time
N = length(idx_start:idx_end);
U_r_SE = squeeze(sqrt(sum(radius_U_SE(:, idx_start:idx_end) .^ 2, 2)))' / N; % RMS of uncertainties

% Create shaded region for uncertainty
curve1 = U_r + U_r_SE; % Upper bound
curve2 = U_r - U_r_SE; % Lower bound
rad2 = [rad, fliplr(rad)]; % X-values for fill
inBetween = [curve1, fliplr(curve2)]; % Y-values for fill

% Plot time-averaged Volume Rate vs. radius
figure("Visible", "off");
fill(rad2, inBetween, Color_std, 'EdgeColor', 'none'); % Shaded region
hold on;
plot(rad, curve1, "Color", Color_std, 'LineWidth', 2); % Upper bound
plot(rad, curve2, "Color", Color_std, 'LineWidth', 2); % Lower bound
plot(rad, U_r, '-k', 'LineWidth', 2); % Mean curve
yline(mean(U_r), '--k', 'LineWidth', 2, 'label', sprintf("%0.2f µL/min", mean(U_r))); % Horizontal mean line

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

ylabel(y_label);
xlabel('Radius (pixels)');
title("Time-Averaged Volume Rate");

% Export plot
exportgraphics(gca, fullfile(ToolBox.path_png, sprintf("%s_mean_%s_radius.png", ToolBox.folder_name, vessel_name)));

% Plot radial variations of Volume Rate over time
figure("Visible", "off");
hold on;

for circleIdx = 1:numCircles
    plot(fullTime, radius_U(circleIdx, :));
end

plot(fullTime, mean(radius_U, 1), 'k', 'LineWidth', 2);

axis padded;
axP = axis;
axis tight;
axT = axis;
axis([axT(1), axT(2), 0, 1.07 * axP(4)]);
box on;
ylabel(y_label);
xlabel('Time (s)');
title("Radial Variations of Volume Rate");
set(gca, 'PlotBoxAspectRatio', [2.5 1 1], 'LineWidth', 2);

% Export plot
exportgraphics(gca, fullfile(ToolBox.path_png, sprintf("%s_variance_%s_time.png", ToolBox.folder_name, vessel_name)));

% Compute total Volume Rate over circles
U_t = squeeze(mean(radius_U, 1)); % Mean over circles
N = size(radius_U_SE, 1);
U_t_SE = squeeze(sqrt(sum(radius_U_SE .^ 2, 1))) / N; % RMS of uncertainties

% Compute statistics for the time range
mean_U = mean(U_t(idx_start:idx_end)); % Time-averaged mean
[max_U, amax] = max(U_t(idx_start:idx_end)); % Maximum value in the range
N = length(idx_start:idx_end);
mean_U_SE = sqrt(sum(U_t_SE(idx_start:idx_end) .^ 2)) / N; % RMS of the uncertainty
[min_U, amin] = min(U_t(idx_start:idx_end)); % Minimum value in the range
max_U_SE = U_t_SE(amax);
min_U_SE = U_t_SE(amin);

% Plot total Volume Rate over time
figure("Visible", "off");
curve1 = U_t + U_t_SE; % Upper bound
curve2 = U_t - U_t_SE; % Lower bound
ft2 = [fullTime, fliplr(fullTime)]; % X-values for fill
inBetween = [curve1, fliplr(curve2)]; % Y-values for fill

hold on;
fill(ft2, inBetween, Color_std, 'EdgeColor', 'none'); % Shaded region
yline(0, 'k-', 'LineWidth', 2); % Zero line
plot(fullTime, curve1, "Color", Color_std, 'LineWidth', 2); % Upper bound
plot(fullTime, curve2, "Color", Color_std, 'LineWidth', 2); % Lower bound
plot(fullTime, U_t, '-k', 'LineWidth', 2); % Mean curve
yline(mean_U, '--k', 'LineWidth', 2); % Horizontal mean line

% Mark the time range used for averaging
plot(fullTime(idx_start), 1.07 * max_U, 'k|', 'MarkerSize', 10, 'LineWidth', 2);
plot(fullTime(idx_end), 1.07 * max_U, 'k|', 'MarkerSize', 10, 'LineWidth', 2);
plot(fullTime(idx_start:idx_end), repmat(1.07 * max_U, 1, idx_end - idx_start + 1), '-k', 'LineWidth', 2);

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

ylabel(y_label);
xlabel('Time (s)');
title(sprintf("Total %s (Avg. %0.2f µL/min)", y_label, mean_U));

% Export plot
exportgraphics(gca, fullfile(ToolBox.path_png, sprintf("%s_allrad_%s_time.png", ToolBox.folder_name, vessel_name)));
exportgraphics(gca, fullfile(ToolBox.path_eps, sprintf("%s_allrad_%s_time.eps", ToolBox.folder_name, vessel_name)));

% Write results to a text file
fileID = fopen(fullfile(ToolBox.path_txt, strcat(ToolBox.folder_name, '_', 'EF_main_outputs', '.txt')), 'a');
fprintf(fileID, 'Flow Rate %s : %f (µL/min) \r\n', vessel_name, mean_U);
fprintf(fileID, 'Flow Rate Standard Deviation %s : %f (µL/min) \r\n', vessel_name, mean_U_SE);
fclose(fileID);

% Write results to json
ToolBox.outputs.(sprintf('FlowRate%s', vessel_name)) = mean_U;
ToolBox.outputs.(sprintf('FlowRateStd%s', vessel_name)) = mean_U_SE;

% New
if contains(vessel_name, 'vein')
    ToolBox.Outputs.add('VenousMeanVolumeRate', mean_U, 'µL/min', mean_U_SE);
    ToolBox.Outputs.add('VenousMaximumVolumeRate', max_U, 'µL/min', max_U_SE);
    ToolBox.Outputs.add('VenousMinimumVolumeRate', min_U, 'µL/min', min_U_SE);
elseif contains(vessel_name, 'artery')
    ToolBox.Outputs.add('ArterialMeanVolumeRate', mean_U, 'µL/min', mean_U_SE);
    ToolBox.Outputs.add('ArterialMinimumVolumeRate', min_U, 'µL/min', min_U_SE);
    ToolBox.Outputs.add('ArterialMaximumVolumeRate', max_U, 'µL/min', max_U_SE);
end

end
