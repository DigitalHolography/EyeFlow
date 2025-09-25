paths = readlines("C:\Users\Mikhalkino\Desktop\in.txt");

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

for ind = 1:length(paths)
    path = paths(ind);

    if isfolder(path)
        path = strcat(path, '\');
    end

    ExecClass = ExecutionClass(path);

    ExecClass.ToolBoxMaster = ToolBoxClass(ExecClass.directory, ExecClass.param_name, 0); %no overwrite

    ExecClass.preprocessData();

    ExecClass.flag_segmentation = 1;
    ExecClass.flag_spectral_analysis = 0;
    ExecClass.flag_bloodFlowVelocity_analysis = 1;
    ExecClass.flag_crossSection_analysis = 1;
    ExecClass.flag_crossSection_figures = 1;

    ExecClass.analyzeData([]);
end

%% Show

Show_multiple_outputs;
