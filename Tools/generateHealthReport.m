function generateHealthReport()

ToolBox = getGlobalToolBox;

% Define file paths
dataFilePath = fullfile(ToolBox.path_txt, strcat(ToolBox.folder_name, '_EF_main_outputs.txt'));
pdfPath = fullfile(ToolBox.path_pdf, sprintf("%s_EyeFlowReport.pdf", ToolBox.folder_name));
directoryName = ToolBox.folder_name;

% Set A4 dimensions in centimeters
a4Width = 21.0; % A4 width in cm
a4Height = 29.7; % A4 height in cm

% Define margins (2 cm on all sides)
margin = 2; % in cm

% Create a new figure for the report
fig = figure('Units', 'centimeters', 'Position', [0 0 a4Width a4Height], ...
    'Color', 'white', 'Visible', 'off', ...
    'PaperSize', [a4Width a4Height], 'PaperPosition', [0 0 a4Width a4Height], ...
    'PaperPositionMode', 'manual'); % Set paper size and position

% Read and parse the data file
data = parseDataFile(dataFilePath);

% Add logo to the top-right corner
addLogo(fig, 'eyeflow_logo.png', a4Width, a4Height);

% Add title with directory name
reportTitle = sprintf('EyeFlow Report - %s', directoryName);
addTitle(fig, reportTitle, margin, a4Width, a4Height);

% Add subtitle with current date/time
currentDateTime = datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss');
subtitleText = sprintf('Generated on: %s', currentDateTime);
addSubtitle(fig, subtitleText, margin, a4Width, a4Height);

% Add data fields
t_factor = ToolBox.stride / ToolBox.fs / 1000;
yPos = 0.85 - (margin / a4Height); % Starting vertical position for annotations (normalized units)
yPos = addField(fig, sprintf('Heart Beat: %.1f bpm', data.heartBeat), yPos, margin, a4Width);
yPos = addField(fig, sprintf('Systoles: %s', mat2str(round(data.systoleIndices * ToolBox.stride / ToolBox.fs / 1000, 2))), yPos, margin, a4Width); % 2 decimal places
yPos = addField(fig, sprintf('Number of Cycles: %d', data.numCycles), yPos, margin, a4Width);
yPos = addField(fig, sprintf('Max Systole Indices: %s', mat2str(round(data.maxSystoleIndices * t_factor, 2))), yPos, margin, a4Width); % 2 decimal places
yPos = addField(fig, sprintf('Min Systole Indices: %s', mat2str(round(data.minSystoleIndices * t_factor, 2))), yPos, margin, a4Width); % 2 decimal places
yPos = addField(fig, sprintf('Mean Flow RateArtery): %.1f ± %.1f µL/min', data.meanBloodVolumeRateArtery, data.stdBloodVolumeRateArtery / 2), yPos, margin, a4Width);
yPos = addField(fig, sprintf('Mean Flow RateVein): %.1f ± %.1f µL/min', data.meanBloodVolumeRateVein, data.stdBloodVolumeRateVein / 2), yPos, margin, a4Width);
yPos = addField(fig, sprintf('Max Systole Flow Rate (Artery): %.1f µL/min', data.maxSystoleBloodVolumeRateArtery), yPos, margin, a4Width);
yPos = addField(fig, sprintf('Min Diastole Flow Rate (Artery): %.1f µL/min', data.minDiastoleBloodVolumeRateArtery), yPos, margin, a4Width);
yPos = addField(fig, sprintf('Stroke Volume (Artery): %.1f nL', data.strokeVolumeArtery), yPos, margin, a4Width);
yPos = addField(fig, sprintf('Total Volume (Artery): %.1f nL', data.totalVolumeArtery), yPos, margin, a4Width);

addImage(fig, fullfile(ToolBox.path_png, 'mask', sprintf("%s_vesselMap.png", directoryName)), yPos, margin, a4Width);

% Save the figure as a PDF using print
print(fig, pdfPath, '-dpdf', '-fillpage'); % Export to PDF with A4 size and margins

% Close the figure
close(fig);
end

% Helper function to parse the data file
function data = parseDataFile(dataFilePath)
% Read the data from the text file
fileContent = fileread(dataFilePath);

% Extract values using regular expressions
data.heartBeat = str2double(regexp(fileContent, 'Heart beat: ([\d.]+)', 'tokens', 'once'));

% Extract Systole Indices and remove trailing comma
systoleIndicesStr = regexp(fileContent, 'Systole Indices: \[([\d,]+)\]', 'tokens', 'once');
data.systoleIndices = str2num(systoleIndicesStr{1}); %#ok<ST2NM>

data.numCycles = str2double(regexp(fileContent, 'Number of Cycles: ([\d]+)', 'tokens', 'once'));

% Extract Max Systole Indices and remove trailing comma
maxSystoleIndicesStr = regexp(fileContent, 'Max Systole Indices: \[([\d,]+)\]', 'tokens', 'once');
data.maxSystoleIndices = str2num(maxSystoleIndicesStr{1}); %#ok<ST2NM>

% Extract Min Systole Indices and remove trailing comma
minSystoleIndicesStr = regexp(fileContent, 'Min Systole Indices: \[([\d,]+)\]', 'tokens', 'once');
data.minSystoleIndices = str2num(minSystoleIndicesStr{1}); %#ok<ST2NM>

data.meanBloodVolumeRateArtery = str2double(regexp(fileContent, 'Flow Rate Artery : ([\d.]+)', 'tokens', 'once'));
data.stdBloodVolumeRateArtery = str2double(regexp(fileContent, 'Flow Rate Standard Deviation Artery : ([\d.]+)', 'tokens', 'once'));
data.meanBloodVolumeRateVein = str2double(regexp(fileContent, 'Flow Rate Vein : ([\d.]+)', 'tokens', 'once'));
data.stdBloodVolumeRateVein = str2double(regexp(fileContent, 'Flow Rate Standard Deviation Vein : ([\d.]+)', 'tokens', 'once'));
data.maxSystoleBloodVolumeRateArtery = str2double(regexp(fileContent, 'MaxSystole Flow Rate Artery : ([\d.]+)', 'tokens', 'once'));
data.minDiastoleBloodVolumeRateArtery = str2double(regexp(fileContent, 'MinDiastole Flow Rate Artery : ([\d.]+)', 'tokens', 'once'));
data.strokeVolumeArtery = str2double(regexp(fileContent, 'Stroke Volume Artery : ([\d.]+)', 'tokens', 'once'));
data.totalVolumeArtery = str2double(regexp(fileContent, 'Total Volume Artery : ([\d.]+)', 'tokens', 'once'));
end

% Helper function to add a title
function addTitle(fig, titleText, margin, a4Width, a4Height)
% Calculate normalized position for the title
titleX = margin / a4Width; % Normalized x position (2 cm margin)
titleY = 1 - (margin / a4Height); % Normalized y position (2 cm margin)

annotation(fig, 'textbox', [titleX titleY - 0.1 1 - 2 * titleX 0.1], 'String', titleText, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
    'FontSize', 20, 'FontWeight', 'bold', 'EdgeColor', 'none', 'Interpreter', 'none');
end

% Helper function to add a data field
function yPos = addField(fig, fieldText, yPos, margin, a4Width)
% Calculate normalized position for the field
fieldX = margin / a4Width; % Normalized x position (2 cm margin)
fieldY = yPos; % Normalized y position (2 cm margin)

annotation(fig, 'textbox', [fieldX fieldY 1 - 2 * fieldX 0.05], 'String', fieldText, ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
    'FontSize', 12, 'EdgeColor', 'none');
yPos = yPos - 0.02; % Move down for the next field
end

% Helper function to add the logo
function addLogo(fig, logoPath, a4Width, a4Height)
% Load the logo image
logoImg = imread(logoPath);

% Calculate normalized position for the logo (top-right corner)
logoWidth = 3; % Width of the logo in cm
logoHeight = size(logoImg, 1) / size(logoImg, 2) * logoWidth; % Maintain aspect ratio
logoX = 1 - (logoWidth / a4Width); % Normalized x position
logoY = 1 - (logoHeight / a4Height); % Normalized y position

% Create axes for the logo
axes('Parent', fig, 'Units', 'normalized', 'Position', [logoX logoY logoWidth / a4Width logoHeight / a4Height]);
h = imshow(logoImg);
axis off;

% Set the transparency (alpha) of the logo to 50%
alpha(h, 0.5);
end

function addSubtitle(fig, subtitleText, margin, a4Width, a4Height)
% Calculate normalized position for the subtitle
subtitleX = margin / a4Width; % Normalized x position (2 cm margin)
subtitleY = 1 - (margin / a4Height) - 0.1; % Normalized y position (below the title)

% Add subtitle with italic font and dark gray color
annotation(fig, 'textbox', [subtitleX subtitleY 1 - 2 * subtitleX 0.05], 'String', subtitleText, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
    'FontSize', 12, 'FontWeight', 'normal', 'FontAngle', 'italic', ...
    'Color', [0.3 0.3 0.3], 'EdgeColor', 'none', 'Interpreter', 'none'); % Dark gray color
end

function yPos = addImage(fig, imagePath, yPos, margin, a4Width)
% ADDIMAGETOPDF Adds an image to the PDF report at specified position
%   Inputs:
%       fig - Figure handle of the report
%       imagePath - Path to the PNG image file
%       position - [x y width height] in centimeters from top-left corner
%       a4Width - Width of A4 page in cm
%       a4Height - Height of A4 page in cm
%
%   Example:
%       addImageToPDF(fig, 'path/to/image.png', [5 10 8 6], a4Width, a4Height)
%       would place the image starting 5cm from left, 10cm from top,
%       with width 8cm and height 6cm

% Check if image file exists
if ~exist(imagePath, 'file')
    error('Image file not found: %s', imagePath);
end

fieldX = margin / a4Width; % Normalized x position (2 cm margin)
fieldY = yPos; % Normalized y position (2 cm margin)

% Create axes for the image
axes('Parent', fig, 'Units', 'normalized', 'Position', [fieldX fieldY 1 - 2 * fieldX 0.05]);

% Read and display the image
try
    img = imread(imagePath);
    imshow(img);
    axis off;

    % Add a subtle border around the image
    rectangle('Position', [0.5, 0.5, size(img, 2) - 1, size(img, 1) - 1], ...
        'EdgeColor', [0.7 0.7 0.7], 'LineWidth', 0.5);
catch ME
    warning(ME.identifier, 'Could not add image to report: %s', ME.message);
    delete(gca); % Remove the axes if image loading failed
end

yPos = yPos - 0.02; % Move down for the next field

end
