clc;
clear;
close all;

%% ---- Step 1: Select and Open TXT File with Folder Paths ----
[txtFile, txtPath] = uigetfile('*.txt', 'Select text file with folder paths');

if isequal(txtFile,0)
    disp('User canceled file selection');
    return;
end

fileID = fopen(fullfile(txtPath, txtFile), 'r');
folderPaths = textscan(fileID, '%s');
fclose(fileID);

folderPaths = folderPaths{1};   % Convert to cell array of strings

%% ---- Step 2: Loop Through Each Folder and Extract Signal ----
figure;
hold on;
grid on;

legendEntries = {};

color1 = [0 0.4470 0.7410];   % MATLAB default blue
color2 = [0.8500 0.3250 0.0980]; % MATLAB default orange

colors = zeros(100, 3);

% First 7 signals -> color1
colors(1:7, :) = repmat(color1, 7, 1);

% Next 12 signals -> color2
colors(8:19, :) = repmat(color2, 12, 1);

% Rest -> color1
colors(20:end, :) = repmat(color1, 100-19, 1);

for i = 1:length(folderPaths)
    
    currentFolder = folderPaths{i};

    [~,foldername,~] = fileparts(currentFolder);
    
    % Define signal file name (CHANGE if needed)

    dir_path_raw = fullfile(currentFolder,"raw");
    h5_files = dir(fullfile(dir_path_raw, '*.h5'));

    RefRawFilePath = fullfile(dir_path_raw, h5_files(1).name);

    % info = h5info(RefRawFilePath, '/');
    % dataset_names = {info.Datasets.Name};

    % M0 = squeeze(h5read(RefRawFilePath, '/moment0'));

    % eyflow_path = fullfile(currentFolder,"eyeflow");

    highest_number = 0; %get_highest_number_in_directories(eyflow_path, "EF");

    eyflow_path = fullfile(currentFolder,"eyeflow",sprintf("%s_EF_%d",foldername,highest_number));

    eyeflow_h5_files = dir(fullfile(eyflow_path,"h5", '*.h5'));

    Refh5FilePath = fullfile(eyflow_path,"h5", eyeflow_h5_files(1).name);

    info = h5info(Refh5FilePath, '/');
    
    signal_M0 = h5read(Refh5FilePath, "/Artery/Velocity/m0SectionSignalRaw/Signal/value");
    
    sysindxlist = h5read(Refh5FilePath, "/Artery/WaveformAnalysis/SystoleIndices/value");
    interpsig = interpSignal(signal_M0, sysindxlist, 512);
    % nnzsection = h5read(Refh5FilePath, "/Artery/Velocity/m0SectionSignalRaw/nnzSection/value");

    data = (interpsig)';

         
    % Check if file has 1 or 2 columns
    if size(data,2) == 1
        y = data;
        x = 1:length(y);
    else
        x = data(:,1);
        y = data(:,2);
    end
    
    plot(x, y, 'LineWidth', 1.2,'Color', colors(i,:));
    
    legendEntries{end+1} = ['Signal ' num2str(i)];
end
    

%% ---- Step 3: Final Plot Adjustments ----
xlabel('Time / Samples');
ylabel('Amplitude');
title('Extracted Signals from All Folders');
legend(legendEntries);
hold off;