function generateA4Report()
% generateA4Report - Creates an A4 PDF with image grid and parameters
%
% Inputs:
%   parameters    - Structure containing parameter names and values

ToolBox = getGlobalToolBox;
path_png = ToolBox.path_png;
path_pdf = ToolBox.path_pdf;
folder_name = ToolBox.folder_name;
outputs = ToolBox.Output;

% Parameters
% Arterious Velocities
parameters.Max_Arterial_Velocity = outputs.ArteryVelocityMax;
parameters.Average_Arterial_Velocity = outputs.ArteryVelocityMean;
parameters.Min_Arterial_Velocity = outputs.ArteryVelocityMin;

% Venous Velocities
parameters.Max_Venous_Velocity = outputs.VeinVelocityMax;
parameters.Average_Venous_Velocity = outputs.VeinVelocityMean;
parameters.Min_Venous_Velocity = outputs.VeinVelocityMin;

% Indexes
parameters.ARI = outputs.ArteryResistivityIndexVelocity;
parameters.VRI = outputs.VeinResistivityIndexVelocity;
parameters.API = outputs.ArteryPulsatilityIndexVelocity;
parameters.VPI = outputs.VeinPulsatilityIndexVelocity;

% Other Parameters
parameters.TimePeakToDescent = outputs.ArteryTimePeakToDescent;
parameters.VeinTimeToPeakFromMin = outputs.VeinTimeToPeakFromMin;
parameters.DicroticNotchVisibility = outputs.DicroticNotchVisibility;

% Volume Rates
parameters.Average_Arterial_Volume_Rate = outputs.ArteryFlowRateMean;
parameters.Average_Venous_Volume_Rate = outputs.VeinFlowRateMean;

% Heart Rate and Blood Pressure
parameters.heart_beat = outputs.HeartBeat;

if isfield(outputs, 'ArterialSystolicFraction')
    parameters.arterial_systolic_fraction = outputs.ArterySystolicFraction;
    parameters.arterial_diastolic_fraction = outputs.ArteryDiastolicFraction;
end

% parameters.time_2_systolic_peak = outputs.Time2SystolicPeak;
parameters.SystoleDuration = outputs.SystoleDuration;
parameters.DiastoleDuration = outputs.DiastoleDuration;
parameters.UnixTimestampFirst = outputs.UnixTimestampFirst;
parameters.UnixTimestampLast = outputs.UnixTimestampLast;

% Create a new figure with A4 paper size (in centimeters)
fig = figure('Units', 'centimeters', ...
    'Position', [0 0 21.0 29.7], ...
    'PaperSize', [21.0 29.7], ...
    'PaperPositionMode', 'auto', ...
    'Color', 'w', ...
    'Visible', 'off');

% === Title ===
t = tiledlayout(fig, 12, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
title(t, folder_name, 'FontSize', 16, 'FontWeight', 'bold', 'Interpreter', 'none');

% Margins
t.Units = 'centimeters';
t.OuterPosition = [1.5 3 18 25.7]; % [left bottom width height]

% === IMAGE GRID ===
vesselTypes = {'artery', 'vein'};

% Create top row (2 columns at 2x height)
for col = 1:2
    name = vesselTypes{col};

    % Top row (1 column at 2x height)
    ax1 = nexttile(t, col, [4 1]);
    vr_combined_path = fullfile(path_png, sprintf('%s_vr_%s.png', folder_name, name));
    vesselmap_path = fullfile(path_png, sprintf('%s_vessel_map_%s.png', folder_name, name));
    mask_path = fullfile(path_png, 'mask', sprintf('%s_M0_%s.png', folder_name, name));
    im1 = loadOrPlaceholder(vr_combined_path, vesselmap_path, mask_path);

    imshow(im1, [], 'Parent', ax1);
    axis(ax1, 'off');

    % Middle row (1 column at 1x height)
    ax2 = nexttile(t, 8 + col, [3 1]);
    vr_path = fullfile(path_png, sprintf('%s_plot_vr_%s.png', folder_name, name));
    ri_path = fullfile(path_png, sprintf('%s_RI_v_%s.png', folder_name, name));
    im2 = loadOrPlaceholder(vr_path, ri_path);

    imshow(im2, [], 'Parent', ax2);
    axis(ax2, 'off');

    % Bottom row (1 column at 1x height)
    ax3 = nexttile(t, 14 + col, [3 1]);
    volume_path = fullfile(path_png, sprintf('%s_strokeAndTotalVolume_%s.png', folder_name, name));

    if strcmp(name, 'artery')
        wave_path = fullfile(path_png, sprintf('%s_ArterialWaveformAnalysis_v_%s.png', folder_name, name));
    else
        wave_path = fullfile(path_png, sprintf('%s_VenousWaveformAnalysis_v_%s.png', folder_name, name));
    end

    find_systole_path = fullfile(path_png, sprintf('%s_find_systoles_indices_%s.png', folder_name, name));

    im3 = loadOrPlaceholder(volume_path, wave_path, find_systole_path);

    imshow(im3, [], 'Parent', ax3);
    axis(ax3, 'off');

end

% === PARAMETERS SECTION ===
axParams = nexttile(t, [2 2]); % last row spans both columns
axis(axParams, 'off');

paramNames = fieldnames(parameters);
numParams = length(paramNames);
paramText = cell(numParams, 1);

for i = 1:numParams
    paramValue = parameters.(paramNames{i});
    fmt = chooseFormat(paramNames{i});
    paramText{i} = sprintf('%s %s', sprintf(fmt, paramValue.value), paramValue.unit);
end

paramFontSize = 12;
paramTitleFontSize = 14;

text(axParams, 0, 1.1, 'Computed Parameters:', ...
    'FontWeight', 'bold', 'FontSize', paramTitleFontSize, ...
    'VerticalAlignment', 'top', 'Interpreter', 'none');

if numParams > 6
    col1 = paramText(1:ceil(numParams / 2));
    col2 = paramText(ceil(numParams / 2) + 1:end);
    text(axParams, 0.0, 0.9, col1, 'VerticalAlignment', 'top', 'FontSize', paramFontSize, 'Interpreter', 'none');
    text(axParams, 0.5, 0.9, col2, 'VerticalAlignment', 'top', 'FontSize', paramFontSize, 'Interpreter', 'none');
else
    text(axParams, 0.0, 0.9, paramText, 'VerticalAlignment', 'top', 'FontSize', paramFontSize, 'Interpreter', 'none');
end

% === Export to PDF ===
set(fig, 'Renderer', 'painters');
output_pdf = fullfile(path_pdf, sprintf('%s_report.pdf', folder_name));
print(fig, output_pdf, '-dpdf', '-r150', '-bestfit');
close(fig);

end

% Helper function to choose the format based on value type
function fmt = chooseFormat(name)

switch name
    case 'UnixTimestampFirst'
        fmt = 'First Unix Timestamp = %d';
    case 'UnixTimestampLast'
        fmt = 'Last Unix Timestamp = %d';
    case 'arterial_systolic_fraction'
        fmt = 'Arterial Systolic Fraction = %.1f';
    case 'arterial_diastolic_fraction'
        fmt = 'Arterial Diastolic Fraction = %.1f';
    case 'time_2_systolic_peak'
        fmt = 'Time to Systolic Peak = %.1f';
    case 'heart_beat'
        fmt = 'HR = %.1f';
    case 'Average_Arterial_Velocity'
        fmt = 'Avg Arterial Velocity = %.2f';
    case 'Max_Arterial_Velocity'
        fmt = 'Max Arterial Velocity = %.2f';
    case 'Max_Venous_Velocity'
        fmt = 'Max Venous Velocity = %.2f';
    case 'Min_Arterial_Velocity'
        fmt = 'Min Arterial Velocity = %.2f';
    case 'Min_Venous_Velocity'
        fmt = 'Min Venous Velocity = %.2f';
    case 'Average_Venous_Velocity'
        fmt = 'Avg Venous Velocity = %.2f';
    case 'TimePeakToDescent'
        fmt = 'Time Peak to Descent = %.2f';
    case 'VeinTimeToPeakFromMin'
        fmt = 'Time to Peak from Min Vein = %.2f';
    case 'DicroticNotchVisibility'
        fmt = 'Dicrotic Notch Visibility = %.0f';
    case 'Average_Arterial_Volume_Rate'
        fmt = 'Avg Arterial Volume Rate = %.2f';
    case 'Average_Venous_Volume_Rate'
        fmt = 'Avg Venous Volume Rate = %.2f';
    case 'SystoleDuration'
        fmt = 'Systole Duration = %.2f';
    case 'DiastoleDuration'
        fmt = 'Diastole Duration = %.2f';
    case 'ARI'
        fmt = 'Arterial Resistivity Index = %.2f';
    case 'VRI'
        fmt = 'Venous Resistivity Index = %.2f';
    case 'API'
        fmt = 'Arterial Pulsatility Index = %.2f';
    case 'VPI'
        fmt = 'Venous Pulsatility Index = %.2f';
    otherwise
        fmt = '%g'; % General format for others
end

end

function im = loadOrPlaceholder(varargin)

for k = 1:numel(varargin)

    if isfile(varargin{k})
        im = imread(varargin{k});
        return;
    end

end

im = ones(200, 200, 3); % White placeholder
end
