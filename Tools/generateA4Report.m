function generateA4Report()
% generateA4Report - Creates an A4 PDF with image grid and parameters
%
% Inputs:
%   parameters    - Structure containing parameter names and values

ToolBox = getGlobalToolBox;
path_png = ToolBox.path_png;
path_pdf = ToolBox.path_pdf;
folder_name = ToolBox.folder_name;
outputs = ToolBox.Outputs;

% Parameters
parameters.average_arterial_velocity = outputs.ArterialMeanVelocity;
parameters.average_venous_velocity = outputs.VenousMeanVelocity;
parameters.TimeToMaxIncreaseSystolic = outputs.TimeToMaxIncreaseSystolic;
parameters.TimePeakToDescent = outputs.TimePeakToDescent;
parameters.TimeToPeakFromMinVein = outputs.TimetoPeakFromMinVein;
parameters.average_arterial_volume_rate = outputs.ArterialMeanVolumeRate;
parameters.average_venous_volume_rate = outputs.VenousMeanVolumeRate;
parameters.heart_beat = outputs.HeartBeat;
parameters.arterial_systolic_fraction = outputs.ArterialSystolicFraction;
parameters.arterial_diastolic_fraction = outputs.ArterialDiastolicFraction;
% parameters.time_2_systolic_peak = outputs.Time2SystolicPeak;
parameters.SystoleDuration = outputs.SystoleDuration;
parameters.DiastoleDuration = outputs.DiastoleDuration;
parameters.ARI = outputs.ArterialResistivityIndexVelocity;
parameters.VRI = outputs.VenousResistivityIndexVelocity;
parameters.API = outputs.ArterialPulsatilityIndexVelocity;
parameters.VPI = outputs.VenousPulsatilityIndexVelocity;

% Create a new figure with A4 paper size (in centimeters)
fig = figure('Units', 'centimeters', 'Position', [0 0 21.0 29.7], ...
    'PaperSize', [21.0 29.7], 'PaperPositionMode', 'auto');

% Set margins (in normalized units)
% topMargin = 0.08;
bottomMargin = 0.10;
leftMargin = 0.05;
% rightMargin = 0.05;

% Add title at the top
axes('Position', [0 0 1 1], 'Visible', 'off');
text(0.5, 1, folder_name, ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'top', ...
    'FontSize', 14, ...
    'FontWeight', 'bold', ...
    'Interpreter', 'none');

% Calculate available space for grid (below title and above parameters)
gridTop = 0.95; % Position below title
gridBottom = bottomMargin + 0.1; % Extra space above parameters
gridHeight = gridTop - gridBottom;

% Define row heights (sum should be 1)
% Top row will be 2x, middle and bottom 1x each (total 4 units)
rowHeights = [0.55, 0.225, 0.21]; % Normalized heights

% Create top row (2 columns at 2x height)
for col = 1:2

    if col == 1
        name = 'artery';
    else
        name = 'vein';
    end

    posX = leftMargin + (col - 1) * 0.45;
    posY = gridBottom + (rowHeights(2) + rowHeights(3)) * gridHeight; % Above middle and bottom rows

    ax = axes('Position', [posX posY 0.45 rowHeights(1) * gridHeight]);
    vr_combined_path = fullfile(path_png, sprintf('%s_combined_vr_%s.png', folder_name, name));

    if isfile(vr_combined_path)
        vr_combined_im = imread(vr_combined_path); % Taller placeholder (2x height)
        imshow(vr_combined_im, []);
    else
        v_path = fullfile(path_png, 'mask', sprintf('%s_vessel_map_%s.png', folder_name, name));
        v_im = imread(v_path); % Taller placeholder (2x height)
        imshow(v_im, []);
    end

    set(ax, 'XTick', [], 'YTick', []);

    % Create middle row (2 columns at 1x height)

    posX = leftMargin + (col - 1) * 0.45;
    posY = gridBottom + rowHeights(3) * gridHeight; % Above bottom row

    ax = axes('Position', [posX posY 0.45 rowHeights(2) * gridHeight]);
    ri_path = fullfile(path_png, sprintf('%s_RI_v_%s.png', folder_name, name));

    if isfile(ri_path)
        ri_im = imread(ri_path); % Standard placeholder
        imshow(ri_im, []);
    else
        placeholder_im = ones(200, 200, 3); % Black placeholder
        imshow(placeholder_im, []);
    end

    set(ax, 'XTick', [], 'YTick', []);

    % Create bottom row (2 columns at 1x height)

    posX = leftMargin + (col - 1) * 0.45;
    posY = gridBottom;

    ax = axes('Position', [posX posY 0.45 rowHeights(3) * gridHeight]);
    volume_path = fullfile(path_png, sprintf('%s_strokeAndTotalVolume_%s.png', folder_name, name));

    if isfile(volume_path)
        volume_im = imread(volume_path); % Standard placeholder
        imshow(volume_im, []);
    else
        if strcmp(name,'artery')
            a_wave_path = fullfile(path_png, sprintf('%s_ArterialWaveformAnalysis_v_%s.png', folder_name, name));
        else
            a_wave_path = fullfile(path_png, sprintf('%s_VenousWaveformAnalysis_v_%s.png', folder_name, name));
        end
        a_wave_im = imread(a_wave_path); % Taller placeholder (2x height)
        imshow(a_wave_im, []);    end
    set(ax, 'XTick', [], 'YTick', []);
end

% Add parameter section at the bottom with more spacing
axes('Position', [leftMargin 0.04 0.8 bottomMargin], 'Visible', 'off');

% Prepare parameter text
paramNames = fieldnames(parameters);
numParams = length(paramNames);
paramText = cell(numParams, 1);

for i = 1:numParams
    paramValue = parameters.(paramNames{i});
    % Format numbers nicely
    if isnumeric(paramValue.value)
        paramText{i} = sprintf('%s: %.2f %s', paramNames{i}, paramValue.value, paramValue.unit);
    else
        paramText{i} = sprintf('%s: %s %s', paramNames{i}, paramValue.value, paramValue.unit);
    end

end

% Split into two columns if more than 6 parameters

% Set larger font size for parameters
paramFontSize = 14; % Increase as desired
paramTitleFontSize = 16;

if numParams > 6
    col1 = paramText(1:ceil(numParams / 2));
    col2 = paramText(ceil(numParams / 2) + 1:end);

    text(0, 1, col1, 'VerticalAlignment', 'top', 'FontSize', paramFontSize, 'Interpreter', 'none');
    text(0.5, 1, col2, 'VerticalAlignment', 'top', 'FontSize', paramFontSize, 'Interpreter', 'none');
else
    text(0, 1, paramText, 'VerticalAlignment', 'top', 'FontSize', paramFontSize, 'Interpreter', 'none');
end

% Add title to parameters section
text(0, 1.15, 'Computed Parameters:', 'FontWeight', 'bold', 'FontSize', paramTitleFontSize);

% Save to PDF
print(fig, fullfile(path_pdf, sprintf('%s_report.pdf', folder_name)), '-dpdf', '-bestfit');

% Close the figure
close(fig);

end
