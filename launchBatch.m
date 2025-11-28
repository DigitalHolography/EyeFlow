function launchBatch()
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

    defaultJson = fullfile(appRoot, "Parameters", "DefaultEyeFlowParamsBatch.json");

    if ~isfile(defaultJson)
        error("CRITICAL: Could not find DefaultEyeFlowParamsBatch.json at: %s", defaultJson);
    else
        fprintf("Parameter file found: %s\n", defaultJson);
    end

    [txt_name, txt_path] = uigetfile('*.txt', 'Select the list of HoloDoppler processed folders');

    if isequal(txt_name, 0)
        fprintf('No file selected. Exiting.\n');
        return;
    end

    fullInputPath = fullfile(txt_path, txt_name);
    paths = strtrim(readlines(fullInputPath));
    paths = paths(paths ~= ""); % remove empty lines

    for ind = 1:length(paths)
        targetPath = fullfile(paths(ind), 'eyeflow');
        jsonDir = fullfile(targetPath, 'json');

        if ~isfolder(jsonDir)
            mkdir(jsonDir);
        end

        % Remove old json files
        delete(fullfile(jsonDir, '*.json'));

        % Copy default parameter file
        copyfile(defaultJson, fullfile(jsonDir, 'input_EF_params.json'));

        % Handle Masks
        maskDir = fullfile(targetPath, 'mask');

        if isfile(fullfile(maskDir, 'forceMaskArtery.png'))
            movefile(fullfile(maskDir, 'forceMaskArtery.png'), fullfile(maskDir, 'oldForceMaskArtery.png'));
        end

        if isfile(fullfile(maskDir, 'forceMaskVein.png'))
            movefile(fullfile(maskDir, 'forceMaskVein.png'), fullfile(maskDir, 'oldForceMaskVein.png'));
        end

    end

    AIModels = AINetworksClass();

    for ind = 1:length(paths)

        p = paths(ind);
        % Ensure path ends with separator
        if isfolder(p) && ~endsWith(p, filesep)
            p = strcat(p, filesep);
        end

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
        MEdisp(e, path);
    end

    try
        ExecClass.Reporter.getA4Report(ME);
        ExecClass.Reporter.saveOutputs();
        ExecClass.Reporter.displayFinalSummary(totalTime);
    catch e
        MEdisp(e, path);
    end

end
