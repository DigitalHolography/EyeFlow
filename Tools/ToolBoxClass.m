classdef ToolBoxClass < handle
% ToolBoxClass holds useful variables for eyeflow processing.

properties
    % Paths
    EF_path char %'X:\...\XXX_HD_X\'
    EF_name char %'XXX_HD_X_EF'
    path_main char %'X:\...\XXX_HD_X\eyeflow'
    path_dir char %'X:\...\XXX_HD_X\eyeflow\XXX_HD_X_EF_X'
    path_png char %'X:\...\XXX_HD_X\eyeflow\XXX_HD_X_EF_X\png'
    path_eps char
    path_pdf char
    path_gif char
    path_txt char
    path_avi char
    path_mp4 char
    path_json char
    path_h5 char
    path_log char
    main_foldername char %'XXX_HD_X'
    param_name char %'input_EF_params.json'
    folder_name char % 'XXX_HD_X_EF_X'
    % Parameters
    stride double
    fs double
    f1 double
    f2 double
    record_time_stamps_us % structure with first and time_stamps
    holo_frames % structure with first and last frames

    % Results
    % Ref % Ref handle to the Execution Class to have access to its properties easily
    Cache % Cache class handle Cache small variables through the execution
    Output % Output class handle Stores outputs through the execution

    params

    params
end

methods

    function obj = ToolBoxClass(path, EF_param_name)
        % Constructor for ToolBoxClass: Initializes paths, parameters, and calculates scaling factors.

        % Store paths and parameters
        obj.EF_path = path;
        obj.param_name = EF_param_name;
        obj.main_foldername = obj.extractFolderName(path);

        % Initialize EyeFlow-related paths
        obj.initializePaths();

        % Load parameters from cache or fall back to defaults
        obj.loadParameters(path);

        % Set up logging (diary)
        obj.setupLogging();

        % Copy input parameters to result folder
        obj.copyInputParameters();

        obj.setGlobalToolBox();
    end

    function mainFolder = extractFolderName(~, path)
        % Helper function to extract the folder name
        split_path = strsplit(path, filesep);
        mainFolder = split_path{end - 1};
    end

    function setGlobalToolBox(obj)
        global ToolBoxGlobal
        ToolBoxGlobal = obj;
    end

    function initializePaths(obj)
        % Helper function to initialize paths for storing eyeflow-related data

        % Define main and subdirectories for storing data
        obj.path_main = fullfile(obj.EF_path, 'eyeflow');
        foldername_EF = strcat(obj.main_foldername, '_EF');

        % Create or identify a unique folder for the current run
        idx = obj.getUniqueFolderIndex(foldername_EF);

        % Special case: if version.txt contains "dev"
        if isfile('version.txt')
            vers = readlines('version.txt');

            if any(contains(vers, 'dev'))
                idx = 0; % Use index 0 for development versions
            end

        end

        % Set the folder name and paths for various data types
        obj.EF_name = foldername_EF;
        obj.folder_name = sprintf('%s_%d', foldername_EF, idx);
        obj.path_dir = fullfile(obj.path_main, obj.folder_name);
        obj.path_png = fullfile(obj.path_dir, 'png');
        obj.path_eps = fullfile(obj.path_dir, 'eps');
        obj.path_pdf = fullfile(obj.path_dir, 'pdf');
        obj.path_txt = fullfile(obj.path_dir, 'txt');
        obj.path_avi = fullfile(obj.path_dir, 'avi');
        obj.path_gif = fullfile(obj.path_dir, 'gif');
        obj.path_mp4 = fullfile(obj.path_dir, 'mp4');
        obj.path_json = fullfile(obj.path_dir, 'json');
        obj.path_h5 = fullfile(obj.path_dir, 'h5');
        obj.path_log = fullfile(obj.path_dir, 'log');

        % Create directories if they don't exist
        obj.createDirectories();
    end

    function idx = getUniqueFolderIndex(obj, folderBaseName)
        % Helper function to determine the unique folder index based on existing directories

        idx = 1;
        list_dir = dir(obj.path_main);

        for i = 1:length(list_dir)

            if contains(list_dir(i).name, folderBaseName)
                match = regexp(list_dir(i).name, '\d+$', 'match');

                if ~isempty(match) && str2double(match{1}) >= idx

                    idx = str2double(match{1}) + 1; % Use the next index

                end

            end

        end

    end

    function createDirectories(obj)
        % Helper function to create necessary directories if they don't exist

        dirs = {obj.path_dir, obj.path_png, obj.path_eps, obj.path_gif, ...
                    obj.path_txt, obj.path_avi, obj.path_mp4, obj.path_json, ...
                    obj.path_log, obj.path_pdf, obj.path_h5};

        for i = 1:length(dirs)

            if ~isfolder(dirs{i})
                mkdir(dirs{i});
            end

        end

    end

    function loadParameters(obj, path)
        % Load or fall back to default parameters from cache or config files

        % Try loading parameters from existing .mat or .json files
        if ~isempty(dir(fullfile(path, ['*', 'RenderingParameters', '*']))) % since HD 2.0
            disp('Reading cache parameters from .json');
            fpath = fullfile(path, dir(fullfile(path, ['*', 'RenderingParameters', '*'])).name);
            decoded_data = jsondecode(fileread(fpath));
            obj.stride = decoded_data.batch_stride;
            obj.fs = decoded_data.fs; % Convert kHz to kHz
            obj.f1 = decoded_data.time_range(1);
            obj.f2 = decoded_data.time_range(2);
            disp('Done.');

        elseif ~isempty(dir(fullfile(path, ['*', 'input_HD_params', '*']))) % since HD 2.9
            disp('Reading cache parameters from input_HD_params');
            fpath = fullfile(path, dir(fullfile(path, ['*', 'input_HD_params', '*'])).name);
            decoded_data = jsondecode(fileread(fpath));
            obj.stride = decoded_data.batch_stride;
            obj.fs = decoded_data.fs; % Convert kHz to kHz
            obj.f1 = decoded_data.time_range(1);

            if isfield(decoded_data, 'record_time_stamps_us')
                obj.record_time_stamps_us = decoded_data.record_time_stamps_us;
            end

            if isfield(decoded_data, 'num_frames')
                obj.holo_frames.first = decoded_data.first_frame;
                obj.holo_frames.last = decoded_data.end_frame;
            end

        elseif isfile(fullfile(path, 'mat', [obj.main_foldername, '.mat']))
            disp('Reading cache parameters from .mat');
            load(fullfile(path, 'mat', [obj.main_foldername, '.mat']), 'cache');
            obj.stride = cache.batch_stride;
            obj.fs = cache.Fs / 1000; % Convert Hz to kHz
            obj.f1 = cache.time_transform.f1;
            obj.f2 = cache.time_transform.f2;

        elseif isfile(fullfile(path, 'Holovibes_rendering_parameters.json'))
            json_txt = fileread(fullfile(path, 'Holovibes_rendering_parameters.json'));
            footer_parsed = jsondecode(json_txt);
            obj.stride = footer_parsed.compute_settings.image_rendering.time_transformation_stride;
            obj.fs = footer_parsed.info.camera_fps / 1000; % Convert FPS to kHz
            obj.f1 = footer_parsed.compute_settings.view.z.start / footer_parsed.compute_settings.image_rendering.time_transformation_size * obj.fs;
            obj.f2 = obj.fs / 2;

        else
            % Default values if no parameters are found
            disp('WARNING: No rendering parameters file found. Using default values.');
            obj.stride = 512; % Default value
            obj.fs = 36; % Default value in kHz
            obj.f1 = 6;
            obj.f2 = obj.fs / 2;
        end

    end

    function setupLogging(obj)
        % Set up logging (diary) for the current session

        diary off % Turn off logging first to avoid logging this script
        diary_filename = fullfile(obj.path_log, sprintf('%s_log.txt', obj.folder_name));
        set(0, 'DiaryFile', diary_filename);
        diary on % Turn on logging
        fprintf("==================================\n");
        fprintf("Current Folder Path: %s\n", obj.EF_path);
        fprintf("Current File: %s\n", obj.folder_name);
        fprintf("Start Time: %s\n", datetime('now', 'Format', 'yyyy/MM/dd HH:mm:ss'));
        fprintf("==================================\n");
    end

    function copyInputParameters(obj)
        % Copy the input parameters to the result folder
        path_dir_json = fullfile(obj.EF_path, 'eyeflow', 'json');
        path_file_json_params_src = fullfile(path_dir_json, obj.param_name);
        path_file_json_params_dest = fullfile(obj.path_json, sprintf("%s_Input_EF_Params.json", obj.folder_name));

        if ~isfile(path_file_json_params_src)
            warning('Source parameter file does not exist: %s', path_file_json_params_src);
        else

            try
                copyfile(path_file_json_params_src, path_file_json_params_dest);
            catch ME
                warning('Failed to copy parameter file: %s\nError: %s', path_file_json_params_src, ME.message);
            end

        end

        % Add the version .txt file in the result folder
        % Resolve path relative to Tools/ToolBoxClass.m -> Parent -> version.txt
        [currentPath, ~, ~] = fileparts(mfilename('fullpath')); 
        [appRoot, ~, ~] = fileparts(currentPath);
        version_file_src = fullfile(appRoot, 'version.txt');
        
        version_file_dest = fullfile(obj.path_dir, sprintf("%s_version.txt", obj.folder_name));
        
        if isfile(version_file_src)
            copyfile(version_file_src, version_file_dest);
        else
            warning('version.txt not found at %s', version_file_src);
        end

        % Add the git .txt file in the result folder
        obj.saveGit();

    end

    function Params = getParams(obj)
        if isempty(obj.params)
            obj.params = Parameters_json(obj.EF_path, obj.param_name);
        end
        Params = obj.params;
    end
    function saveGit(obj)
        % SAVING GIT VERSION
        % In the txt file in the folder : "log"

        % Get the current branch name
        gitBranchCommand = 'git symbolic-ref --short HEAD';
        [statusBranch, resultBranch] = system(gitBranchCommand);

        if statusBranch == 0
            resultBranch = strtrim(resultBranch);
            MessBranch = 'Current branch : %s \r';
        else
            vers = readlines('version.txt');
            MessBranch = sprintf('PulseWave GitHub version %s\r', char(vers));
        end

        % Get the latest commit hash
        gitHashCommand = 'git rev-parse HEAD';
        [statusHash, resultHash] = system(gitHashCommand);

        if statusHash == 0 %hash command was successful
            resultHash = strtrim(resultHash);
            MessHash = 'Latest Commit Hash : %s \r';
        else
            MessHash = '';
        end

        % Get the most recent tag
        gitTagCommand = 'git describe --tags';
        [statusTag, resultTag] = system(gitTagCommand);

        if statusTag == 0 %tag command was successful
            resultTag = strtrim(resultTag);
            MessTag = 'Most recent tag : %s \r';
        else
            MessTag = '';
        end

        % Print the results
        fprintf('----------------------------------\rGIT VERSION :\r');
        fprintf(MessBranch, resultBranch);
        fprintf(MessHash, resultHash);
        fprintf(MessTag, resultTag);
        fprintf('----------------------------------\r');

        % Save the results to a text file
        logFilePath = fullfile(obj.path_dir, sprintf('%s_git_version.txt', obj.folder_name));
        fileID = fopen(logFilePath, 'w');

        if fileID == -1
            error('Cannot open file: %s', logFilePath);
        end

        fprintf(fileID, '----------------------------------\rGIT VERSION :\r');
        fprintf(fileID, MessBranch, resultBranch);
        fprintf(fileID, MessHash, resultHash);
        fprintf(fileID, MessTag, resultTag);
        fprintf(fileID, '----------------------------------\r');
        fclose(fileID);

    end

    function ToolBox = setCache(obj, Cache)
        obj.Cache = Cache;
        ToolBox = obj;
    end

    function ToolBox = setOutput(obj, Output)
        obj.Output = Output;
        ToolBox = obj;
    end

end

end
