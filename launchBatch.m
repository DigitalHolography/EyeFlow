%% Get the list of holodoppler folders to be processed from user
fprintf("=== EYEFLOWPROCESS START ===\n");

beginComputerTime = sprintf("Eyeprocess Begin Computer Time: %s\n", datetime('now', 'Format', 'yyyy/MM/dd HH:mm:ss'));

version_tag = readlines("version.txt");

fprintf("EyeFlow version : %s\n", version_tag);
fprintf("\n%s\n",beginComputerTime);

[txt_name, txt_path] = uigetfile('*.txt', 'Select the list of HoloDoppler processed folders');

if isequal(txt_name, 0)
    fprintf('No file selected. Exiting.\n');
    return;
end

paths = strtrim(readlines(fullfile(txt_path, txt_name)));
paths = paths(paths ~= ""); % remove empty lines

% Add necessary paths
addpath("BloodFlowVelocity\", "BloodFlowVelocity\Elastography\", "CrossSection\", ...
    "Loading\", "Parameters\", "Preprocessing\", ...
    "Scripts\", "Segmentation\", "SHAnalysis\", "Tools\", "Outputs\");
%% ensure set default parameters and no forced mask

for ind = 1:length(paths)
    path = fullfile(paths(ind), 'eyeflow');

    if ~isfolder(fullfile(path, 'json'))
        mkdir(fullfile(path, 'json'));
    end

    delete(fullfile(fullfile(path, 'json'), '*.json')); % remove old json files
    copyfile(fullfile('Parameters', 'DefaultEyeFlowParams.json'), fullfile(path, 'json', 'input_EF_params.json'));

    if isfile(fullfile(path, 'mask', 'forceMaskArtery.png'))
        movefile(fullfile(path, 'mask', 'forceMaskArtery.png'), fullfile(path, 'mask', 'oldForceMaskArtery.png'));
    end

    if isfile(fullfile(path, 'mask', 'forceMaskVein.png'))
        movefile(fullfile(path, 'mask', 'forceMaskVein.png'), fullfile(path, 'mask', 'oldForceMaskVein.png'));
    end

end

%% launch

% Generate timestamped log file name
t = datetime('now', 'Format', 'yyyyMMdd_HHmmss');
logFileName = sprintf('log_%s.txt', char(t));

if ~isfolder("Logs")
    mkdir("Logs");
end

fprintf('Log saving to Eyeflow\%s\n', fullfile('Logs', logFileName));
fid = fopen(fullfile('Logs', logFileName), 'a'); % 'a' for append if needed

AIModels = AINetworksClass();

for ind = 1:length(paths)

    path = paths(ind);
    fprintf(fid, 'Execution of Eyeflow routine on %s  ;  %d/%d\n', path, ind, length(path));

    
    

    if isfolder(path)
        path = strcat(path, '\');
    end

    tic;
    runAnalysisBlock(path, AIModels);
    ti = toc;
    fprintf(fid, 'Execution time: %.2f seconds\n\n', ti);
end

fclose(fid);

endComputerTime = sprintf("Eyeprocess End Computer Time: %s\n", datetime('now', 'Format', 'yyyy/MM/dd HH:mm:ss'));


fprintf('Log saved to Eyeflow\%s\n', fullfile('Logs', logFileName));

fprintf("\n   (. ❛ ᴗ ❛.)\n");
fprintf("\n");

fprintf("\n%s\n",beginComputerTime);
fprintf("\n%s\n",endComputerTime);
fprintf("\n");

fprintf("=== EYEFLOWPROCESS END ===\n");

% %% Show
% try
%     ShowOutputs(paths, 'Logs');
% catch ME
%     MEdisp(ME, 'Logs');
% end

%%

function runAnalysisBlock(path, AIModels)

totalTime = tic;

ME = [];

try
    ExecClass = ExecutionClass(path);

    ExecClass.AINetworks = AIModels;

    ExecClass.ToolBoxMaster = ToolBoxClass(ExecClass.directory, ExecClass.param_name);

    ExecClass.preprocessData();

    ExecClass.flag_segmentation = 0;
    ExecClass.flag_spectral_analysis = 0;
    ExecClass.flag_bloodFlowVelocity_analysis = 0;
    ExecClass.flag_pulseWaveVelocity = 0;
    ExecClass.flag_crossSection_analysis = 0;
    ExecClass.flag_crossSection_export = 0;

    ExecClass.analyzeData([]);
catch e
    ME = e;
    MEdisp(e, path);
end


ReporterTimer = tic;
fprintf("\n----------------------------------\n" + ...
    "Generating Reports\n" + ...
"----------------------------------\n");
ExecClass.Reporter.getA4Report(ME);
ExecClass.Reporter.saveOutputs();
fprintf("- Reporting took : %ds\n", round(toc(ReporterTimer)))
ExecClass.Reporter.displayFinalSummary(totalTime);

end
