paths = readlines("Y:\todo\Test_list\segmentation_test_list.txt");
% Add necessary paths
addpath("BloodFlowVelocity\", "BloodFlowVelocity\Elastography\", "CrossSection\", ...
    "Loading\", "Parameters\", "Preprocessing\", ...
    "Scripts\", "Segmentation\", "SHAnalysis\", "Tools\");
%% ensure set default parameters and no forced mask

% for ind = 1:length(paths)
%     split_path = strsplit(paths(ind), '\');
%     main_foldername = split_path{end};
%     folder_name = strcat(main_foldername, '_EF');
%     path = fullfile(paths(ind),'eyeflow');
%     list_dir = dir(path);
%     idx = 0;
%     for i=1:length(list_dir)
%         if contains(list_dir(i).name, folder_name)
%             match = regexp(list_dir(i).name, '\d+$', 'match');
%             if ~isempty(match) && str2double(match{1}) >= idx
%                 idx = str2double(match{1}); %suffix
%             end
%         end
%     end
%     last_folder_name = sprintf('%s_%d', folder_name, idx);
%
%     copyfile(fullfile('Parameters','DefaultEyeFlowParams.json'),fullfile(path,'json','input_EF_params.json'));
%
%     if isfile(fullfile(path,'mask','forceMaskArtery.png'))
%         movefile(fullfile(path,'mask','forceMaskArtery.png'),fullfile(path,'mask','oldForceMaskArtery.png'));
%     end
% end

%% launch

% Generate timestamped log file name
t = datetime('now', 'Format', 'yyyyMMdd_HHmmss');
logFileName = sprintf('log_%s.txt', char(t));

if ~isfolder("Logs")
    mkdir("Logs");
end

fprintf('Log saving to %s :\n', logFileName);
fid = fopen(fullfile('Logs', logFileName), 'a'); % 'a' for append if needed

for ind = 1:length(paths)

    path = paths(ind);
    fprintf(fid, 'Execution of Eyeflow routine on %s  ;  %d/%d\n', path, ind, length(path));

    if isfolder(path)
        path = strcat(path, '\');
    end

    tic;
    runAnalysisBlock(path);
    ti = toc;
    fprintf(fid, 'Execution time: %.2f seconds\n\n', ti);
end
fclose(fid);

fprintf('Log saved to %s\n', logFileName);

%% Show
try
    ShowOutputs(paths, 'Logs');
catch ME
    MEdisp(ME, 'Logs')
end

%%

function runAnalysisBlock(path)

try
    ExecClass = ExecutionClass(path);
    ExecClass.ToolBoxMaster = ToolBoxClass(ExecClass.directory, ExecClass.param_name, 0); %no overwrite

    ExecClass.preprocessData();

    ExecClass.flag_segmentation = 1;
    ExecClass.flag_spectral_analysis = 0;
    ExecClass.flag_bloodFlowVelocity_analysis = 0;
    ExecClass.flag_crossSection_analysis = 0;
    ExecClass.flag_crossSection_figures = 0;

    ExecClass.analyzeData([]);
catch ME
    MEdisp(ME, path)
end

end
