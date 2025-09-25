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
parameters.Max_Arterial_Velocity = outputs.ArterialMaximumVelocity;
parameters.Average_Arterial_Velocity = outputs.ArterialMeanVelocity;
parameters.Min_Arterial_Velocity = outputs.ArterialMinimumVelocity;

% Venous Velocities
parameters.Max_Venous_Velocity = outputs.VenousMaximumVelocity;
parameters.Average_Venous_Velocity = outputs.VenousMeanVelocity;
parameters.Min_Venous_Velocity = outputs.VenousMinimumVelocity;

% Indexes
parameters.ARI = outputs.ArterialResistivityIndexVelocity;
parameters.VRI = outputs.VenousResistivityIndexVelocity;
parameters.API = outputs.ArterialPulsatilityIndexVelocity;
parameters.VPI = outputs.VenousPulsatilityIndexVelocity;

% Other Parameters
parameters.TimePeakToDescent = outputs.TimePeakToDescent;
parameters.TimeToPeakFromMinVein = outputs.TimetoPeakFromMinVein;
parameters.DicroticNotchVisibility = outputs.DicroticNotchVisibility;

% Volume Rates
parameters.Average_Arterial_Volume_Rate = outputs.ArterialMeanVolumeRate;
parameters.Average_Venous_Volume_Rate = outputs.VenousMeanVolumeRate;

% Heart Rate and Blood Pressure
parameters.heart_beat = outputs.HeartBeat;

parameters.arterial_systolic_fraction = outputs.ArterialSystolicFraction;
parameters.arterial_diastolic_fraction = outputs.ArterialDiastolicFraction;
% parameters.time_2_systolic_peak = outputs.Time2SystolicPeak;
parameters.SystoleDuration = outputs.SystoleDuration;
parameters.DiastoleDuration = outputs.DiastoleDuration;
parameters.UnixTimestampFirst = outputs.UnixTimestampFirst;
parameters.UnixTimestampLast = outputs.UnixTimestampLast;

% Create a new figure with A4 paper size (in centimeters)
fig = figure('Units', 'centimeters', 'Position', [0 0 21.0 29.7], ...
    'PaperSize', [21.0 29.7], 'PaperPositionMode', 'auto');

% Set margins (in normalized units)
% topMargin = 0.08;
bottomMargin = 0.20;
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
rowHeights = [0.50, 0.25, 0.25]; % Normalized heights

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

        if strcmp(name, 'artery')
            a_wave_path = fullfile(path_png, sprintf('%s_ArterialWaveformAnalysis_v_%s.png', folder_name, name));
        else
            a_wave_path = fullfile(path_png, sprintf('%s_VenousWaveformAnalysis_v_%s.png', folder_name, name));
        end

        a_wave_im = imread(a_wave_path); % Taller placeholder (2x height)
        imshow(a_wave_im, []);
    end

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
    fmt = chooseFormat(paramNames{i});

    % Format numbers nicely
    paramText{i} = sprintf('%s %s', sprintf(fmt, paramValue.value), paramValue.unit);

end

% Split into two columns if more than 6 parameters

% Set larger font size for parameters
paramFontSize = 12; % Increase as desired
paramTitleFontSize = 14;

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

% Set figure renderer to control compression
set(fig, 'Renderer', 'painters'); % Vector graphics where possible

% Save to PDF with controlled resolution
print(fig, fullfile(path_pdf, sprintf('%s_report.pdf', folder_name)), ...
    '-dpdf', '-r150', '-bestfit', '-image'); % -r150 sets 150 DPI

% Close the figure
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
    case 'TimeToPeakFromMinVein'
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
