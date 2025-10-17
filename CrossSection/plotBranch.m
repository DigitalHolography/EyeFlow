function plotBranch(branch_U, branch_U_SE, fullTime, idx_start, idx_end, vessel_name, type_name)

% Get global ToolBox and parameters
ToolBox = getGlobalToolBox;

numBranches = size(branch_U, 1); % Number of branches
bId = 1:numBranches; % Branch indices
N = length(idx_start:idx_end);

% Define color for shaded regions
Color_std = [0.7 0.7 0.7];

if strcmp(type_name, 'velocity')
    title_str = 'Velocity';
    y_label = 'Velocity (mm/s)';
else
    title_str = 'flow_rate';
    y_label = 'Flow Rate (µL/min)';
end

% Compute time-averaged mean and standard deviation
U_b = squeeze(mean(branch_U(:, idx_start:idx_end), 2))'; % Mean over time
U_b_SE = squeeze(sqrt(sum(branch_U_SE(:, idx_start:idx_end) .^ 2, 2)))' / N; % RMS of uncertainties

% Plot time-averaged vs. branch
figure("Visible", "off");
hold on
errorbar(bId, U_b, U_b_SE, 'k', "LineStyle", "none", 'LineWidth', 2); % Error bars
bar(bId, U_b, 1, "FaceColor", Color_std, 'EdgeColor', 'k', 'LineWidth', 2); % Bar plot

% Format plot
axis padded;
axP = axis;
axis tight;
axT = axis;
axis([axT(1), axT(2), 0, 1.07 * axP(4)]);

box on
set(gca, 'LineWidth', 2);
set(gca, 'PlotBoxAspectRatio', [1.618, 1, 1])
ax = gca;
ax.LineStyleOrderIndex = 1; % Reset if needed
ax.SortMethod = 'depth'; % Try changing sorting method
ax.Layer = 'top'; % This may help in some cases

ylabel(y_label);
xlabel('branch id.');

% Export plot
exportgraphics(gca, fullfile(ToolBox.path_png, sprintf("%s_branch_errorbar_%s_%s.png", ToolBox.folder_name, title_str, vessel_name)));

% Plot branch variations of Flow Rate over time
figure("Visible", "off");
hold on;

for branchIdx = 1:numBranches
    plot(fullTime, branch_U(branchIdx, :));
end

plot(fullTime, mean(branch_U, 1), 'k', 'LineWidth', 2);

axis padded;
axP = axis;
axis tight;
axT = axis;
axis([axT(1), axT(2), 0, 1.07 * axP(4)]);
box on;
ylabel(y_label);
xlabel('Time (s)');
set(gca, 'PlotBoxAspectRatio', [2.5 1 1], 'LineWidth', 2);

% Export plot
exportgraphics(gca, fullfile(ToolBox.path_png, sprintf("%s_branch_%s_variance_%s.png", ToolBox.folder_name, title_str, vessel_name)));

end
