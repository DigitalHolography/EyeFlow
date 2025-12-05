function launchBatch(varargin)
    % Initialize Application
    appRoot = fileparts(mfilename('fullpath'));
    versionFile = fullfile(appRoot, "version.txt");

    if isfile(versionFile)
        version_tag = readlines(versionFile);
        fprintf("EyeFlow version : %s\n", version_tag);
    else
        fprintf("EyeFlow version : Unknown (version.txt not found at %s)\n", versionFile);
    end

    beginComputerTime = datetime('now', 'Format', 'yyyy/MM/dd HH:mm:ss');
    fprintf("Begin Computer Time: %s\n", beginComputerTime);

    if ~isdeployed
        % In Development: Add all subfolders to path
        fprintf("Running in Development Mode\n");
        addpath(fullfile(appRoot, "BloodFlowVelocity"));
        addpath(fullfile(appRoot, "BloodFlowVelocity", "Elastography"));
        addpath(fullfile(appRoot, "CrossSection"));
        addpath(fullfile(appRoot, "Loading"));
        addpath(fullfile(appRoot, "Parameters"));
        addpath(fullfile(appRoot, "Preprocessing"));
        addpath(fullfile(appRoot, "Scripts"));
        addpath(fullfile(appRoot, "Segmentation"));
        addpath(fullfile(appRoot, "SHAnalysis"));
        addpath(fullfile(appRoot, "Tools"));
        addpath(fullfile(appRoot, "Outputs"));
    else
        fprintf("Running in Deployed Mode.\n");
    end

    defaultJson = fullfile(appRoot, "Parameters", "DefaultEyeFlowParamsBatch.json");

    if ~isfile(defaultJson)
        error("CRITICAL: Could not find DefaultEyeFlowParamsBatch.json at: %s", defaultJson);
    else
        fprintf("Parameter file found: %s\n", defaultJson);
    end

    % --- ARGUMENT PARSING & PATH SELECTION ---
    paths = string.empty;
    mode = 'single'; % default to single holo mode
    inputPath = "";

    % Check arguments
    if nargin > 0
        arg1 = string(varargin{1});

        if arg1 == "-b"
            mode = 'batch';

            if nargin > 1
                inputPath = string(varargin{2});
            end

        else
            mode = 'single';
            inputPath = arg1;
        end

    end

    try

        if strcmp(mode, 'batch')
            % --- BATCH MODE ---
            if inputPath == ""
                [txt_name, txt_path] = uigetfile('*.txt', 'Select the list of HoloDoppler processed folders');

                if isequal(txt_name, 0)
                    fprintf('No file selected. Exiting.\n');
                    return;
                end

                fullInputPath = fullfile(txt_path, txt_name);
            else
                fullInputPath = inputPath;

                if ~isfile(fullInputPath)
                    error("Batch file not found: %s", fullInputPath);
                end

            end

            fprintf("Running in Batch Mode: %s\n", fullInputPath);
            paths = strtrim(readlines(fullInputPath));
            paths = paths(paths ~= ""); % remove empty lines

        else
            % --- SINGLE HOLO MODE ---
            if inputPath == ""
                [holo_name, holo_path] = uigetfile('*.holo', 'Select the .holo file');

                if isequal(holo_name, 0)
                    fprintf('No file selected. Exiting.\n');
                    return;
                end

                fullInputPath = fullfile(holo_path, holo_name);
            else
                fullInputPath = inputPath;

                if ~isfile(fullInputPath)
                    error("Holo file not found: %s", fullInputPath);
                end

            end

            fprintf("Running in Single Mode: %s\n", fullInputPath);

            % Find the latest HD folder for this holo file
            latestHD = findLatestHDFolder(fullInputPath);

            if ismissing(latestHD)
                fprintf("No HoloDoppler (HD) folder found for: %s\n", fullInputPath);
                return;
            end

            fprintf("Selected latest HD folder: %s\n", latestHD);
            paths = [latestHD];
        end

    catch ME
        fprintf("Error during path selection: %s\n", ME.message);
        return;
    end

    if isempty(paths)
        fprintf("No paths to process. Exiting.\n");
        return;
    end

    % --- PRE-PROCESSING (JSON & MASKS) ---
    fprintf("Preparing %d folder(s)...\n", length(paths));

    for ind = 1:length(paths)
        targetPath = paths(ind);

        % check for empty lines if readlines failed to catch them
        if targetPath == "" || ismissing(targetPath); continue; end

        jsonDir = fullfile(targetPath, 'eyeflow', 'json');

        if ~isfolder(jsonDir)
            mkdir(jsonDir);
        end

        % Remove old json files
        delete(fullfile(jsonDir, '*.json'));

        % Copy default parameter file
        copyfile(defaultJson, fullfile(jsonDir, 'input_EF_params.json'));

        % Handle Masks
        maskDir = fullfile(targetPath, 'eyeflow', 'mask');

        efPath = fullfile(targetPath, 'eyeflow');
        jsonDir = fullfile(efPath, 'json');

        if ~isfolder(jsonDir)
            mkdir(jsonDir);
        end

        delete(fullfile(jsonDir, '*.json'));
        copyfile(defaultJson, fullfile(jsonDir, 'input_EF_params.json'));

        maskDir = fullfile(efPath, 'mask');

        if isfile(fullfile(maskDir, 'forceMaskArtery.png'))
            movefile(fullfile(maskDir, 'forceMaskArtery.png'), fullfile(maskDir, 'oldForceMaskArtery.png'));
        end

        if isfile(fullfile(maskDir, 'forceMaskVein.png'))
            movefile(fullfile(maskDir, 'forceMaskVein.png'), fullfile(maskDir, 'oldForceMaskVein.png'));
        end

    end

    % --- AI LOADING ---
    fprintf("Loading AI Models...\n");
    AIModels = AINetworksClass();

    % --- EXECUTION LOOP ---
    for ind = 1:length(paths)
        p = paths(ind);

        if p == "" || ismissing(p); continue; end

        % Ensure path ends with separator
        if isfolder(p) && ~endsWith(p, filesep)
            p = strcat(p, filesep);
        end

        fprintf("Processing: %s\n", p);
        runAnalysisBlock(p, AIModels);
    end

    endComputerTime = datetime('now', 'Format', 'yyyy/MM/dd HH:mm:ss');
    fprintf("All done! Total time elapsed: %s\n", string(endComputerTime - beginComputerTime));
end

function runAnalysisBlock(path, AIModels)
    totalTime = tic;
    ME = [];

    try
        ExecClass = ExecutionClass(path);

        ExecClass.AINetworks = AIModels;

        % Ensure ExecutionClass uses the directory passed to it
        ExecClass.ToolBoxMaster = ToolBoxClass(ExecClass.directory, ExecClass.param_name);

        ExecClass.preprocessData();

        ExecClass.flag_segmentation = 1;
        ExecClass.flag_spectral_analysis = 1;
        ExecClass.flag_bloodFlowVelocity_analysis = 1;
        ExecClass.flag_pulseWaveVelocity = 1;
        ExecClass.flag_crossSection_analysis = 1;
        ExecClass.flag_crossSection_export = 1;

        ExecClass.analyzeData([]);
    catch e
        ME = e;
        % Check if MEdisp exists, otherwise use disp
        if exist('MEdisp', 'file')
            MEdisp(e, path);
        else
            disp(e.message);
        end

    end

    try
        ExecClass.Reporter.getA4Report(ME);
        ExecClass.Reporter.saveOutputs();
        ExecClass.Reporter.displayFinalSummary(totalTime);
    catch e

        if exist('MEdisp', 'file')
            MEdisp(e, path);
        else
            disp(e.message);
        end

    end

end

function latestHDPath = findLatestHDFolder(holoPath)
    % Finds the HD folder with the highest render number in the same directory as the holo file
    % Naming convention: [holo_filename]_HD_[RENDER_NUMBER]

    latestHDPath = string(missing);

    [parentDir, fname, ~] = fileparts(holoPath);

    % Get list of folders in parent directory
    contents = dir(parentDir);
    dirFlags = [contents.isdir];
    subFolders = contents(dirFlags);

    maxRenderNum = -1;

    % Escaped pattern for regex (escape special regex characters in filename)
    escapedName = regexptranslate('escape', fname);
    pattern = "^" + escapedName + "_HD_(\d+)$";

    for k = 1:length(subFolders)
        thisName = string(subFolders(k).name);

        tokens = regexp(thisName, pattern, 'tokens');

        if ~isempty(tokens)
            % Extract the number
            renderNum = str2double(tokens{1}{1});

            if renderNum > maxRenderNum
                maxRenderNum = renderNum;
                latestHDPath = fullfile(parentDir, thisName);
            end

        end

    end

end
