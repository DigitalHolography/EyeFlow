classdef eyeflow < matlab.apps.AppBase

% Properties that correspond to app components – now PRIVATE
properties (Access = private)
    EyeFlowUIFigure matlab.ui.Figure
    RootGrid matlab.ui.container.GridLayout

    ReferenceDirectory matlab.ui.control.TextArea
    statusLamp matlab.ui.control.Lamp

    % Top Buttons
    LoadFolderButton matlab.ui.control.Button
    LoadHoloButton matlab.ui.control.Button
    ClearButton matlab.ui.control.Button
    FolderManagementButton matlab.ui.control.Button

    % Third Row Buttons
    EditMasksButton matlab.ui.control.Button
    EditParametersButton matlab.ui.control.Button
    PlayMomentsButton matlab.ui.control.Button

    % Fourth Row Buttons
    OpenDirectoryButton matlab.ui.control.Button
    ReProcessButton matlab.ui.control.Button

    % Checkboxes
    segmentationCheckBox matlab.ui.control.CheckBox
    bloodFlowAnalysisCheckBox matlab.ui.control.CheckBox
    pulseVelocityCheckBox matlab.ui.control.CheckBox
    generateCrossSectionSignalsCheckBox matlab.ui.control.CheckBox
    exportCrossSectionResultsCheckBox matlab.ui.control.CheckBox
    spectralAnalysisCheckBox matlab.ui.control.CheckBox

    % Execute
    NumberofWorkersSpinner matlab.ui.control.Spinner
    NumberofWorkersSpinnerLabel matlab.ui.control.Label
    ExecuteButton matlab.ui.control.Button
    ImageDisplay matlab.ui.control.Image
end

% Public data (kept as originally)
properties (Access = public)
    file ExecutionClass
    AINetworks
    drawer_list = {}
end

methods (Access = public)

    % Set the ImageDisplay image source with the same pre‑processing as Load()
    function setDisplayImage(app, imageMatrix)
        % imageMatrix – 2D or 3D numeric array
        if ndims(imageMatrix) == 3
            grayImg = mean(imageMatrix, 3);
        else
            grayImg = imageMatrix;
        end

        rgbImg = repmat(rescale(grayImg), [1 1 3]);
        [numX, numY] = size(rgbImg);
        app.ImageDisplay.ImageSource = imresize(rgbImg, [max(numX, numY) max(numX, numY)]);
    end

    % Provide access to the main figure (needed by FolderManagementUI)
    function fig = getMainFigure(app)
        fig = app.EyeFlowUIFigure;
    end

    function value = getWidgetValue(app, widgetName)
        % Return the current value of the UI widget named widgetName.
        w = app.(widgetName);

        if isempty(w)
            value = [];
        elseif isa(w, 'matlab.ui.control.CheckBox')
            value = logical(w.Value);
        elseif isa(w, 'matlab.ui.control.DropDown')
            value = w.Value;
        elseif isa(w, 'matlab.ui.control.ListBox')
            value = w.Value;
        elseif isa(w, 'matlab.ui.control.NumericEditField') || isa(w, 'matlab.ui.control.Spinner')
            value = double(w.Value);
        else
            value = w.Value;
        end

    end

    function setWidgetValue(app, widgetName, value)
        % Set the value of the UI widget named widgetName.
        w = app.(widgetName);
        if isempty(w), return; end

        if isa(w, 'matlab.ui.control.CheckBox')
            w.Value = logical(value);
        elseif isa(w, 'matlab.ui.control.DropDown')

            if ismember(value, w.Items)
                w.Value = value;
            else
                w.Value = w.Items{1};
            end

        elseif isa(w, 'matlab.ui.control.ListBox')
            allItems = w.Items;

            if iscell(value)
                w.Value = intersect(value, allItems, 'stable');
            else
                w.Value = {};
            end

        elseif isa(w, 'matlab.ui.control.NumericEditField') || isa(w, 'matlab.ui.control.Spinner')
            w.Value = double(value(1));
        else
            w.Value = value;
        end

    end

    function Load(app, path)
        % Update lamp color to indicate loading
        app.statusLamp.Color = [1, 1/2, 0]; % Orange
        drawnow;

        if isfolder(path)
            path = strcat(path, '\');
        end

        totalLoadingTime = tic;

        try
            % Add file
            app.file = ExecutionClass(path);

            app.file.AINetworks = app.AINetworks;

            % Compute the mean of M0 along the third dimension
            mean_M0 = mean(app.file.M0, 3);
            % Display the mean image in the uiimage component
            img = repmat(rescale(mean_M0), [1 1 3]);
            [numX, numY] = size(img);
            app.ImageDisplay.ImageSource = imresize(img, [max(numX, numY) max(numX, numY)]); % Rescale the image for display

            % Enable buttons
            app.ExecuteButton.Enable = true;
            app.ClearButton.Enable = true;
            app.EditParametersButton.Enable = true;
            app.EditMasksButton.Enable = true;
            app.PlayMomentsButton.Enable = true;
            app.OpenDirectoryButton.Enable = true;
            app.ReProcessButton.Enable = true;
            app.ReferenceDirectory.Value = path;

            % Update lamp color to indicate success
            app.statusLamp.Color = [0, 1, 0]; % Green

        catch ME
            MEdisp(ME, path);
            diary off
            % Update lamp color to indicate error
            app.statusLamp.Color = [1, 0, 0]; % Red
        end

        % Update checkbox states after loading
        app.CheckboxValueChanged();

        fprintf("----------------------------------\n")
        fprintf("- Total Load timing took : %ds\n", round(toc(totalLoadingTime)))
    end

end

% Callbacks that handle component events
methods (Access = public)

    % Code that executes after component creation
    function startupFcn(app)

        % Close any existing parallel pool
        delete(gcp('nocreate'))

        if exist("version.txt", 'file')
            v = readlines('version.txt');
            fprintf("==================================\n " + ...
                "Welcome to EyeFlow %s\n" + ...
                "----------------------------------\n" + ...
                "Developed by the DigitalHolographyFoundation\n" + ...
                "==================================\n", v(1));
        end

        % Add necessary paths
        addpath("BloodFlowVelocity\", "BloodFlowVelocity\Elastography\", "CrossSection\", ...
            "Loading\", "Parameters\", "Preprocessing\", "Outputs\", ...
            "Scripts\", "Segmentation\", "SHAnalysis\", "Tools\");

        % Set the UI title
        app.EyeFlowUIFigure.Name = ['EyeFlow ', char(v(1))];

        % Display splash screen
        displaySplashScreen();

        % Initialize checkbox states
        app.CheckboxValueChanged();
        set(groot, 'defaultFigureColor', 'w');
        set(groot, 'defaultAxesFontSize', 14);
        set(groot, 'DefaultTextFontSize', 10); % For text objects (e.g., annotations)

        app.AINetworks = AINetworksClass();
    end

    function LoadFromTxt(app)

        [selected_file, path] = uigetfile('*.txt');

        if (selected_file)
            files_lines = readlines(fullfile(path, selected_file));

            for nn = 1:length(files_lines)

                if ~isempty(files_lines(nn))
                    app.drawer_list{end + 1} = files_lines(nn);
                end

            end

        end

    end

    % Button pushed function: LoadfolderButton
    function LoadfolderButtonPushed(app, ~)

        if ~isempty(app.file)
            last_dir = app.file.directory;
        else
            last_dir = [];
        end

        selected_dir = uigetdir(last_dir);

        if selected_dir == 0
            fprintf(2, 'No folder selected\n');
            return
        else
            % Clearing before loading
            ClearButtonPushed(app);
            app.Load(selected_dir);
        end

    end

    % Button pushed function: LoadHoloButton
    function LoadHoloButtonPushed(app, ~)

        [selected_holo, path_holo] = uigetfile('*.holo');

        if selected_holo == 0
            fprintf(2, 'No file selected\n');
        else
            % Clearing before loading
            ClearButtonPushed(app);
            app.Load(fullfile(path_holo, selected_holo));
        end

    end

    % Button pushed function: ExecuteButton
    function err = ExecuteButtonPushed(app, ~)

        err = [];

        if isempty(app.file)
            fprintf(2, "No input loaded.\n")
            return
        end

        % Update lamp color to indicate processing
        app.statusLamp.Color = [1, 1/2, 0]; % Orange
        drawnow;

        % Actualizes the input Parameters

        fprintf("\n==================================\n");

        fclose all; % Close any open files

        app.file.params_names = checkEyeFlowParamsFromJson(app.file.directory); % checks compatibility between found EF params and Default EF params of this version of EF.
        params = Parameters_json(app.file.directory, app.file.params_names{1});

        if params.json.NumberOfWorkers > 0 && params.json.NumberOfWorkers < app.NumberofWorkersSpinner.Limits(2)
            app.NumberofWorkersSpinner.Value = params.json.NumberOfWorkers;
        end

        for i = 1:length(app.file.params_names)

            app.file.param_name = app.file.params_names{i};

            totalTime = tic;

            fprintf("\n==================================\n")

            app.file.flag_segmentation = app.segmentationCheckBox.Value;
            app.file.flag_bloodFlowVelocity_analysis = app.bloodFlowAnalysisCheckBox.Value;
            app.file.flag_pulseWaveVelocity = app.pulseVelocityCheckBox.Value;
            app.file.flag_crossSection_analysis = app.generateCrossSectionSignalsCheckBox.Value;
            app.file.flag_crossSection_export = app.exportCrossSectionResultsCheckBox.Value;
            app.file.flag_spectral_analysis = app.spectralAnalysisCheckBox.Value;

            err = [];

            try
                app.file.ToolBoxMaster = ToolBoxClass(app.file.directory, app.file.param_name);

                if ~app.file.is_preprocessed
                    parfor_arg = app.NumberofWorkersSpinner.Value;
                    setupParpool(parfor_arg);
                    app.file.preprocessData();
                end

                app.file.analyzeData(app);

                % Update lamp color to indicate success
                app.statusLamp.Color = [0, 1, 0]; % Green

            catch ME
                err = ME;
                MEdisp(ME, app.file.directory);

                % Update lamp color to indicate warning
                app.statusLamp.Color = [1, 0, 0]; % Red
                diary off
            end

            % Generate reports and outputs

            ReporterTimer = tic;
            fprintf("\n----------------------------------\n" + ...
                "Generating Reports\n" + ...
            "----------------------------------\n");
            app.file.Reporter.getA4Report(err);
            app.file.Reporter.saveOutputs();
            fprintf("- Reporting took : %ds\n", round(toc(ReporterTimer)))
            app.file.Reporter.displayFinalSummary(totalTime);

        end

        % Update checkbox states after execution
        app.CheckboxValueChanged();

    end

    function PlayMomentsButtonPushed(app, ~)

        try

            if app.file.is_preprocessed
                disp('inputs after preprocess.')
            else
                disp('inputs before preprocess.')
            end

            implay(rescale(app.file.M0));
            implay(rescale(app.file.M1));
            implay(rescale(app.file.M2));
        catch
            disp('Input not well loaded')
        end

    end

    % Button pushed function: ClearButton
    function ClearButtonPushed(app, ~)

        if ~isempty(app.file)
            clear app.file;
        end

        close all;

        app.ReferenceDirectory.Value = "";

        app.ExecuteButton.Enable = false;
        app.ClearButton.Enable = false;
        app.EditParametersButton.Enable = false;
        app.OpenDirectoryButton.Enable = false;
        app.ReProcessButton.Enable = false;
        app.EditMasksButton.Enable = false;
        app.PlayMomentsButton.Enable = false;

        % Update checkbox states
        app.CheckboxValueChanged();
    end

    % Callback function for Open Directory button
    function OpenDirectoryButtonPushed(app, ~)

        try
            % Open the directory in the system file explorer
            winopen(fullfile(app.file.directory, 'eyeflow')); % For Windows
            % Use `open(app.file.directory)` for macOS/Linux
        catch
            fprintf(2, "No valid directory loaded.\n");
        end

    end

    % Callback function for Open Directory button
    function ReProcessButtonPushed(app, ~)

        % Update lamp color to indicate processing
        app.statusLamp.Color = [1, 1/2, 0]; % Orange
        drawnow;

        try
            parfor_arg = app.NumberofWorkersSpinner.Value;
            setupParpool(parfor_arg);
            app.file = ExecutionClass(app.file.directory);
            app.file.AINetworks = app.AINetworks;
            app.file.preprocessData();

            % Update lamp color to indicate success
            app.statusLamp.Color = [0, 1, 0]; % Green
            drawnow;
        catch ME
            MEdisp(ME, app.file.directory);

            % Update lamp color to indicate warning
            app.statusLamp.Color = [1, 0, 0]; % Red
            drawnow;
            diary off
        end

    end

    function FolderManagementButtonPushed(app, ~)
        persistent fmUI

        if isempty(fmUI) || ~isvalid(fmUI)
            fmUI = FolderManagementUI(app);
        else
            fmUI.bringToFront();
        end

    end

    % Button pushed function: EditParametersButton
    function EditParametersButtonPushed(app, ~)

        try
            main_path = fullfile(app.file.directory, 'eyeflow');

            if isfile(fullfile(main_path, 'json', app.file.param_name))
                disp(['opening : ', fullfile(main_path, 'json', app.file.param_name)])
                winopen(fullfile(main_path, 'json', app.file.param_name));
            else
                disp(['couldn''t open : ', fullfile(main_path, 'json', app.file.param_name)])
            end

        catch
            fprintf(2, "no input loaded\n")
        end

    end

    % Button pushed function: EditMasksButton
    function EditMasksButtonPushed(app, ~)
        ToolBox = getGlobalToolBox;

        if isempty(ToolBox) || ~strcmp(app.file.directory, ToolBox.EF_path)
            ToolBox = ToolBoxClass(app.file.directory, app.file.param_name);
        end

        if ~isempty(app.file)

            if ~isfolder(fullfile(ToolBox.EF_path, 'mask'))
                mkdir(fullfile(ToolBox.EF_path, 'mask'))
            end

            try
                winopen(fullfile(ToolBox.EF_path, 'mask'));
            catch
                disp("opening failed.")
            end

            if ~app.file.is_preprocessed
                parfor_arg = app.NumberofWorkersSpinner.Value;
                setupParpool(parfor_arg);
                app.file = app.file.preprocessData();
            end

            try
                list_dir = dir(ToolBox.EF_path);

                for i = 1:length(list_dir)

                    if contains(list_dir(i).name, ToolBox.folder_name)
                        match = regexp(list_dir(i).name, '\d+$', 'match');

                        if ~isempty(match) && str2double(match{1}) >= idx
                            idx = str2double(match{1}); %suffix
                        end

                    end

                end

                path_dir = fullfile(ToolBox.EF_path, ToolBox.folder_name);

                disp(['Copying from : ', fullfile(path_dir, 'png', 'mask')])
                copyfile(fullfile(path_dir, 'png', 'mask', sprintf("%s_maskArtery.png", ToolBox.folder_name)), fullfile(ToolBox.EF_path, 'mask', 'MaskArtery.png'));
                copyfile(fullfile(path_dir, 'png', 'mask', sprintf("%s_maskVein.png", ToolBox.folder_name)), fullfile(ToolBox.EF_path, 'mask', 'MaskVein.png'));
            catch
                disp("last auto mask copying failed.")
            end

            try

                copyfile(fullfile(ToolBox.EF_path, 'png', sprintf("%s_M0.png", ToolBox.folder_name)), fullfile(ToolBox.EF_path, 'mask', 'M0.png'));
                folder_name = ToolBox.folder_name;
                list_dir = dir(ToolBox.EF_path);
                idx = 0;

                for i = 1:length(list_dir)

                    if contains(list_dir(i).name, folder_name)
                        match = regexp(list_dir(i).name, '\d+$', 'match');

                        if ~isempty(match) && str2double(match{1}) >= idx
                            idx = str2double(match{1}); %suffix
                        end

                    end

                end

                folder_name = sprintf('%s_%d', ToolBox.folder_name, idx);
                copyfile(fullfile(ToolBox.EF_path, folder_name, 'gif', sprintf("%s_M0.gif", folder_name)), fullfile(ToolBox.EF_path, 'mask', 'M0.gif'));
            catch

                disp("last M0 png and gif copying failed")
            end

            try

                copyfile(fullfile(ToolBox.EF_path, 'png', sprintf("%s_M0.png", ToolBox.folder_name)), fullfile(ToolBox.EF_path, 'mask', 'M0.png'));
                folder_name = ToolBox.folder_name;
                list_dir = dir(ToolBox.EF_path);
                idx = 0;

                for i = 1:length(list_dir)

                    if contains(list_dir(i).name, folder_name)
                        match = regexp(list_dir(i).name, '\d+$', 'match');

                        if ~isempty(match) && str2double(match{1}) >= idx
                            idx = str2double(match{1}); %suffix
                        end

                    end

                end

                folder_name = sprintf('%s_%d', ToolBox.folder_name, idx);
                copyfile(fullfile(ToolBox.EF_path, folder_name, 'png', 'mask', sprintf("%s_DiaSysRGB.png", folder_name)), fullfile(ToolBox.EF_path, 'mask', 'DiaSysRGB.png'));
            catch

                disp("Diasys png failed")

            end

            % try
            % %   Commented until further fixes MESSAGE TO ZACHARIE
            %     openmaskinpaintnet(fullfile(ToolBox.EF_path,'mask','M0.png'), fullfile(ToolBox.EF_path,'mask','DiaSysRGB.png'));
            % catch
            %     disp("paint.net macro failed")
            % end

        else

            fprintf(2, "no input loaded\n")

        end

    end

    function CheckboxValueChanged(app, ~)
        % Callback function triggered when any checkbox is clicked.
        % This function enforces the rules for enabling/disabling checkboxes.

        % Segmentation Checkbox is always enabled
        app.segmentationCheckBox.Enable = true;

        % Determine the current state of the file analysis
        is_segmented = false;
        is_velocityAnalyzed = false;
        is_volumeRateAnalyzed = false;

        if ~isempty(app.file)
            is_segmented = app.file.is_segmented;
            is_velocityAnalyzed = app.file.is_velocityAnalyzed;
            is_volumeRateAnalyzed = app.file.is_volumeRateAnalyzed;
        end

        % Enable/disable bloodFlowAnalysisCheckBox, pulseVelocityCheckBox, and spectralAnalysisCheckBox
        if app.segmentationCheckBox.Value || is_segmented
            app.bloodFlowAnalysisCheckBox.Enable = true;
            app.spectralAnalysisCheckBox.Enable = true;
        else
            app.bloodFlowAnalysisCheckBox.Enable = false;
            app.spectralAnalysisCheckBox.Enable = false;
            app.bloodFlowAnalysisCheckBox.Value = false; % Turn off if disabled
            app.spectralAnalysisCheckBox.Value = false;
        end

        % Enable/disable generateCrossSectionSignalsCheckBox
        if app.bloodFlowAnalysisCheckBox.Value || is_velocityAnalyzed
            app.generateCrossSectionSignalsCheckBox.Enable = true;
            app.pulseVelocityCheckBox.Enable = true;
        else
            app.pulseVelocityCheckBox.Enable = false;
            app.generateCrossSectionSignalsCheckBox.Enable = false;
            app.pulseVelocityCheckBox.Value = false;
            app.generateCrossSectionSignalsCheckBox.Value = false; % Turn off if disabled
        end

        % Enable/disable exportCrossSectionResultsCheckBox
        if app.generateCrossSectionSignalsCheckBox.Value || is_volumeRateAnalyzed
            app.exportCrossSectionResultsCheckBox.Enable = true;
        else
            app.exportCrossSectionResultsCheckBox.Enable = false;
            app.exportCrossSectionResultsCheckBox.Value = false; % Turn off if disabled
        end

    end

end

% =========================================================================
% Component initialisation – modernised with panels and grid layouts
% =========================================================================
methods (Access = private)

    function createComponents(app)
        % Master entry point – builds all UI components
        pathToMLAPP = fileparts(mfilename('fullpath'));

        app.createFigure(pathToMLAPP);
        app.createRootGrid();
        app.createFileSelectionPanel();
        app.createProcessingOptionsPanel();
        app.createExecutionToolsPanel();
        app.createImageDisplay();

        % Global appearance
        fontname(app.EyeFlowUIFigure, 'Arial');
        fontsize(app.EyeFlowUIFigure, 12, "points");

        app.EyeFlowUIFigure.Visible = 'on';
    end

    % ---------------------------------------------------------------------
    function createFigure(app, pathToMLAPP)
        screenSize = get(0, 'ScreenSize');
        appWidth = 1050;
        appHeight = 500;
        appX = (screenSize(3) - appWidth) / 2;
        appY = (screenSize(4) - appHeight) / 2;

        app.EyeFlowUIFigure = uifigure('Visible', 'off', ...
            'Position', [appX, appY, appWidth, appHeight], ...
            'Color', [0.2 0.2 0.2], ...
            'Name', 'EyeFlow', ...
            'Icon', fullfile(pathToMLAPP, 'eyeflow_logo.png'), ...
            'WindowStyle', 'normal');
    end

    % ---------------------------------------------------------------------
    function createRootGrid(app)
        % Main grid: two columns (left panels, right image)
        app.RootGrid = uigridlayout(app.EyeFlowUIFigure, [3, 2], ...
            'ColumnWidth', {'fit', '1x'}, ...
            'RowHeight', {'1x', 'fit', 'fit'}, ...
            'Padding', [10 10 10 10], ...
            'RowSpacing', 10, ...
            'ColumnSpacing', 10, ...
            'BackgroundColor', [0.2 0.2 0.2]);
    end

    % ---------------------------------------------------------------------
    function createFileSelectionPanel(app)
        import matlab.ui.container.Panel
        import matlab.ui.container.GridLayout
        import matlab.ui.control.*

        backgroundColor = [0.2 0.2 0.2];
        darkBackgroundColor = [0.15 0.15 0.15];
        fontColor = [1 1 1];
        grayButtonColor = [0.5 0.5 0.5];

        panel = uipanel(app.RootGrid, ...
            'Title', 'File Selection', ...
            'BackgroundColor', backgroundColor, ...
            'ForegroundColor', fontColor, ...
            'BorderType', 'line', ...
            'FontWeight', 'bold');
        panel.Layout.Row = 1;
        panel.Layout.Column = 1;

        grid = uigridlayout(panel, [3, 4], ...
            'ColumnWidth', {'fit', 'fit', 'fit', 'fit'}, ...
            'RowHeight', {'fit', '1x', 'fit'}, ...
            'Padding', [10 10 10 10], ...
            'RowSpacing', 5, ...
            'ColumnSpacing', 5, ...
            'BackgroundColor', backgroundColor);

        % Row 1: Load buttons
        app.LoadFolderButton = uibutton(grid, 'push', ...
            'Text', 'Load Folder', ...
            'BackgroundColor', [0.5 0.5 0.5], ...
            'FontColor', [1 1 1], ...
            'ButtonPushedFcn', @(~, ~) app.LoadfolderButtonPushed());
        app.LoadFolderButton.Layout.Row = 1;
        app.LoadFolderButton.Layout.Column = 1;

        app.LoadHoloButton = uibutton(grid, 'push', ...
            'Text', 'Load Holo', ...
            'BackgroundColor', grayButtonColor, ...
            'FontColor', [1 1 1], ...
            'ButtonPushedFcn', @(~, ~) app.LoadHoloButtonPushed());
        app.LoadHoloButton.Layout.Row = 1;
        app.LoadHoloButton.Layout.Column = 2;

        app.ClearButton = uibutton(grid, 'push', ...
            'Text', 'Clear', ...
            'BackgroundColor', grayButtonColor, ...
            'FontColor', [1 1 1], ...
            'Enable', 'off', ...
            'ButtonPushedFcn', @(~, ~) app.ClearButtonPushed());
        app.ClearButton.Layout.Row = 1;
        app.ClearButton.Layout.Column = 3;

        app.FolderManagementButton = uibutton(grid, 'push', ...
            'Text', 'Folder Management', ...
            'BackgroundColor', grayButtonColor, ...
            'FontColor', [1 1 1], ...
            'ButtonPushedFcn', @(~, ~) app.FolderManagementButtonPushed());
        app.FolderManagementButton.Layout.Row = 1;
        app.FolderManagementButton.Layout.Column = 4;

        % Row 2: Directory display + lamp
        app.ReferenceDirectory = uitextarea(grid, ...
            'Value', '', ...
            'Editable', 'off', ...
            'BackgroundColor', darkBackgroundColor, ...
            'FontColor', [1 1 1]);
        app.ReferenceDirectory.Layout.Row = 2;
        app.ReferenceDirectory.Layout.Column = [1 4];

        % Row 3: Edit tools buttons
        app.EditParametersButton = uibutton(grid, 'push', ...
            'Text', 'Edit Parameters', ...
            'BackgroundColor', grayButtonColor, ...
            'FontColor', [1 1 1], ...
            'Enable', 'off', ...
            'Tooltip', 'Open the JSON parameter file', ...
            'ButtonPushedFcn', @(~, ~) app.EditParametersButtonPushed());
        app.EditParametersButton.Layout.Row = 3;
        app.EditParametersButton.Layout.Column = 1;

        app.EditMasksButton = uibutton(grid, 'push', ...
            'Text', 'Edit Masks', ...
            'BackgroundColor', grayButtonColor, ...
            'FontColor', [1 1 1], ...
            'Enable', 'off', ...
            'Tooltip', 'Open mask folder for manual editing', ...
            'ButtonPushedFcn', @(~, ~) app.EditMasksButtonPushed());
        app.EditMasksButton.Layout.Row = 3;
        app.EditMasksButton.Layout.Column = 2;

        app.PlayMomentsButton = uibutton(grid, 'push', ...
            'Text', 'Play Moments', ...
            'BackgroundColor', grayButtonColor, ...
            'FontColor', [1 1 1], ...
            'Enable', 'off', ...
            'Tooltip', 'Play M0/M1/M2 videos', ...
            'ButtonPushedFcn', @(~, ~) app.PlayMomentsButtonPushed());
        app.PlayMomentsButton.Layout.Row = 3;
        app.PlayMomentsButton.Layout.Column = 3;
    end

    % ---------------------------------------------------------------------
    function createProcessingOptionsPanel(app)
        panel = uipanel(app.RootGrid, ...
            'Title', 'Processing Options', ...
            'BackgroundColor', [0.2 0.2 0.2], ...
            'ForegroundColor', [1 1 1], ...
            'BorderType', 'line', ...
            'FontWeight', 'bold');
        panel.Layout.Row = 2;
        panel.Layout.Column = 1;

        grid = uigridlayout(panel, [6, 2], ...
            'ColumnWidth', {'1x', '1x'}, ...
            'RowHeight', repmat({'fit'}, 1, 6), ...
            'Padding', [10 10 10 10], ...
            'RowSpacing', 8, ...
            'BackgroundColor', [0.2 0.2 0.2]);

        % Segmentation (full width)
        app.segmentationCheckBox = uicheckbox(grid, ...
            'Text', 'Segmentation', ...
            'Value', true, ...
            'FontColor', [1 1 1], ...
            'Tooltip', 'Segment vessels using U-Net', ...
            'ValueChangedFcn', @(~, ~) app.CheckboxValueChanged());
        app.segmentationCheckBox.Layout.Row = 1;
        app.segmentationCheckBox.Layout.Column = [1 2];

        % Blood Flow Analysis & Pulse Analysis
        app.bloodFlowAnalysisCheckBox = uicheckbox(grid, ...
            'Text', 'Blood Flow Analysis', ...
            'Value', true, ...
            'FontColor', [1 1 1], ...
            'Tooltip', 'Compute blood flow velocity', ...
            'ValueChangedFcn', @(~, ~) app.CheckboxValueChanged());
        app.bloodFlowAnalysisCheckBox.Layout.Row = 2;
        app.bloodFlowAnalysisCheckBox.Layout.Column = 1;

        app.pulseVelocityCheckBox = uicheckbox(grid, ...
            'Text', 'Pulse Analysis', ...
            'Value', false, ...
            'FontColor', [1 1 1], ...
            'Tooltip', 'Flexural pulse wave velocity', ...
            'ValueChangedFcn', @(~, ~) app.CheckboxValueChanged());
        app.pulseVelocityCheckBox.Layout.Row = 2;
        app.pulseVelocityCheckBox.Layout.Column = 2;

        % Cross Section Signals & Export
        app.generateCrossSectionSignalsCheckBox = uicheckbox(grid, ...
            'Text', 'Generate Cross Section Signals', ...
            'Value', false, ...
            'FontColor', [1 1 1], ...
            'Tooltip', 'Blood flow profiles across vessels', ...
            'ValueChangedFcn', @(~, ~) app.CheckboxValueChanged());
        app.generateCrossSectionSignalsCheckBox.Layout.Row = 3;
        app.generateCrossSectionSignalsCheckBox.Layout.Column = 1;

        app.exportCrossSectionResultsCheckBox = uicheckbox(grid, ...
            'Text', 'Export Cross Section Results', ...
            'Value', false, ...
            'FontColor', [1 1 1], ...
            'Tooltip', 'Export cross-section analysis', ...
            'ValueChangedFcn', @(~, ~) app.CheckboxValueChanged());
        app.exportCrossSectionResultsCheckBox.Layout.Row = 3;
        app.exportCrossSectionResultsCheckBox.Layout.Column = 2;

        % Spectral Analysis (full width)
        app.spectralAnalysisCheckBox = uicheckbox(grid, ...
            'Text', 'Spectral Analysis', ...
            'Value', false, ...
            'FontColor', [1 1 1], ...
            'Tooltip', 'Spectral analysis of blood flow', ...
            'ValueChangedFcn', @(~, ~) app.CheckboxValueChanged());
        app.spectralAnalysisCheckBox.Layout.Row = 4;
        app.spectralAnalysisCheckBox.Layout.Column = [1 2];
    end

    % ---------------------------------------------------------------------
    function createExecutionToolsPanel(app)
        panel = uipanel(app.RootGrid, ...
            'Title', 'Execution & Tools', ...
            'BackgroundColor', [0.2 0.2 0.2], ...
            'ForegroundColor', [1 1 1], ...
            'BorderType', 'line', ...
            'FontWeight', 'bold');
        panel.Layout.Row = 3;
        panel.Layout.Column = 1;

        grid = uigridlayout(panel, [2, 4], ...
            'ColumnWidth', {'1x', 'fit', 'fit', 'fit'}, ...
            'RowHeight', {'fit', 'fit'}, ...
            'Padding', [10 10 10 10], ...
            'RowSpacing', 8, ...
            'ColumnSpacing', 8, ...
            'BackgroundColor', [0.2 0.2 0.2]);

        % Execute button (green)
        app.ExecuteButton = uibutton(grid, 'push', ...
            'Text', 'Execute', ...
            'BackgroundColor', [0.2 0.6 0.2], ...
            'FontColor', [1 1 1], ...
            'Enable', 'off', ...
            'ButtonPushedFcn', @(~, ~) app.ExecuteButtonPushed());
        app.ExecuteButton.Layout.Row = 1;
        app.ExecuteButton.Layout.Column = 1;

        % Open Directory button
        app.OpenDirectoryButton = uibutton(grid, 'push', ...
            'Text', 'Open Directory', ...
            'BackgroundColor', [0.5 0.5 0.5], ...
            'FontColor', [1 1 1], ...
            'Enable', 'off', ...
            'ButtonPushedFcn', @(~, ~) app.OpenDirectoryButtonPushed());
        app.OpenDirectoryButton.Layout.Row = 1;
        app.OpenDirectoryButton.Layout.Column = 2;

        % Preprocess button
        app.ReProcessButton = uibutton(grid, 'push', ...
            'Text', 'Preprocess', ...
            'BackgroundColor', [0.5 0.5 0.5], ...
            'FontColor', [1 1 1], ...
            'Enable', 'off', ...
            'ButtonPushedFcn', @(~, ~) app.ReProcessButtonPushed());
        app.ReProcessButton.Layout.Row = 1;
        app.ReProcessButton.Layout.Column = 3;

        % Workers label and spinner
        app.NumberofWorkersSpinnerLabel = uilabel(grid, ...
            'Text', 'Workers:', ...
            'HorizontalAlignment', 'right', ...
            'FontColor', [1 1 1]);
        app.NumberofWorkersSpinnerLabel.Layout.Row = 1;
        app.NumberofWorkersSpinnerLabel.Layout.Column = 4;

        app.NumberofWorkersSpinner = uispinner(grid, ...
            'Limits', [0 parcluster('local').NumWorkers], ...
            'Value', min(10, floor(parcluster('local').NumWorkers / 2)), ...
            'BackgroundColor', [0.15 0.15 0.15], ...
            'FontColor', [1 1 1]);
        app.NumberofWorkersSpinner.Layout.Row = 2;
        app.NumberofWorkersSpinner.Layout.Column = 4;

        app.statusLamp = uilamp(grid, ...
            'Color', [0 1 0]); % green
        app.statusLamp.Layout.Row = 2;
        app.statusLamp.Layout.Column = 1;
    end

    % ---------------------------------------------------------------------
    function createImageDisplay(app)
        app.ImageDisplay = uiimage(app.RootGrid);
        app.ImageDisplay.Layout.Row = [1 3];
        app.ImageDisplay.Layout.Column = 2;
        app.ImageDisplay.ScaleMethod = 'fit';
        app.ImageDisplay.BackgroundColor = [0.15 0.15 0.15];
    end

end

% App creation and deletion
methods (Access = public)

    % Construct app
    function app = eyeflow

        % Create UIFigure and components
        createComponents(app)

        % Register the app with App Designer
        registerApp(app, app.EyeFlowUIFigure)

        % Execute the startup function
        runStartupFcn(app, @startupFcn)

        if nargout == 0
            clear app
        end

    end

    % Code that executes before app deletion
    function delete(app)

        % Delete UIFigure when app is deleted
        delete(app.EyeFlowUIFigure)
    end

end

end
