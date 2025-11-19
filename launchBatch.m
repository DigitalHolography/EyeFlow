function launchBatch()

fprintf("=== EYEFLOWPROCESS START ===\n");

% In Dev: Returns the folder containing this script (e.g., C:\Projects\EyeFlow)
% In Deployed: Returns the extraction folder (e.g., ...\Cache\EyeFlo2\EyeFlowProcess)
appRoot = fileparts(mfilename('fullpath'));

try
    cd(appRoot);
    fprintf("Current Working Directory changed to: %s\n", pwd);
catch ME
    warning("Could not change directory to appRoot");
end

fprintf("Application Root: %s\n", appRoot);
beginComputerTime = sprintf("Eyeprocess Begin Computer Time: %s\n", datetime('now', 'Format', 'yyyy/MM/dd HH:mm:ss'));

% Locate version.txt relative to the script
versionFile = fullfile(appRoot, "version.txt");
if isfile(versionFile)
    version_tag = readlines(versionFile);
    fprintf("EyeFlow version : %s\n", version_tag);
else
    fprintf("EyeFlow version : Unknown (version.txt not found at %s)\n", versionFile);
end

fprintf("\n%s\n", beginComputerTime);

% Manage paths
if ~isdeployed
    % In Development: Add all subfolders to path
    fprintf("Running in Development Mode - Adding paths...\n");
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

% We expect the 'Parameters' folder to be a sibling of this script in the cache.
defaultJson = fullfile(appRoot, "Parameters", "DefaultEyeFlowParamsBatch.json");

% Fallback check (compiler sometimes shifts structure by one level)
if ~isfile(defaultJson)
    fallbackPath = fullfile(fileparts(appRoot), "Parameters", "DefaultEyeFlowParamsBatch.json");
    if isfile(fallbackPath)
        defaultJson = fallbackPath;
    else
        error("CRITICAL: Could not find DefaultEyeFlowParamsBatch.json at: %s", defaultJson);
    end
end
fprintf("Parameter file found: %s\n", defaultJson);

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

% Save logs in the same folder as the selected text file (Persistent location)
logDir = fullfile(txt_path, 'Logs'); 
if ~isfolder(logDir)
    mkdir(logDir);
end

t = datetime('now', 'Format', 'yyyyMMdd_HHmmss');
logFileName = sprintf('log_%s.txt', char(t));
logFullPath = fullfile(logDir, logFileName);

fprintf("Log saving to: %s\n", logFullPath);
fid = fopen(logFullPath, 'a'); 

if fid == -1
    error("Could not open log file for writing at: %s", logFullPath);
end

AIModels = AINetworksClass();

for ind = 1:length(paths)

    p = paths(ind);
    fprintf(fid, 'Execution of Eyeflow routine on %s  ;  %d/%d\n', p, ind, length(paths));

    % Ensure path ends with separator
    if isfolder(p) && ~endsWith(p, filesep)
        p = strcat(p, filesep);
    end

    tic;
    runAnalysisBlock(p, AIModels);
    ti = toc;
    fprintf(fid, 'Execution time: %.2f seconds\n\n', ti);
end

fclose(fid);

endComputerTime = sprintf("Eyeprocess End Computer Time: %s\n", datetime('now', 'Format', 'yyyy/MM/dd HH:mm:ss'));

fprintf('Log saved to: %s\n', logFullPath);
fprintf("\n   (. ❛ ᴗ ❛.)\n");
fprintf("\n%s\n", beginComputerTime);
fprintf("\n%s\n", endComputerTime);
fprintf("=== EYEFLOWPROCESS END ===\n");

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

ReporterTimer = tic;
fprintf("\n----------------------------------\n" + ...
    "Generating Reports\n" + ...
"----------------------------------\n");

% Pass error stack if any
ExecClass.Reporter.getA4Report(ME);
ExecClass.Reporter.saveOutputs();

fprintf("- Reporting took : %ds\n", round(toc(ReporterTimer)))
ExecClass.Reporter.displayFinalSummary(totalTime);

end