classdef FolderManagementUI < handle
% FolderManagementUI   Batch folder manager for EyeFlow
%
%   FM = FolderManagementUI(APP) creates a non-modal figure that allows
%   managing a list of EyeFlow measurement folders and running batch
%   analysis.

properties (Access = private)
    MainApp % Reference to the main EyeFlow app
    TextArea % uitextarea displaying the folder list
    WorkerSpinner % uispinner for number of parallel workers
    Figure % Allow reading from outside, but only class can set it
end

methods

    function obj = FolderManagementUI(mainApp)

        arguments
            mainApp
        end

        obj.MainApp = mainApp;
        obj.createUI();

        fontname(obj.Figure, 'Arial');
        fontsize(obj.Figure, 12, "points");
    end

    function bringToFront(obj)
        % Bring the existing UI figure to the front
        if ~isempty(obj.Figure) && isvalid(obj.Figure)
            figure(obj.Figure);
        end

    end

end

methods (Access = private)

    function createUI(obj)
        app = obj.MainApp;

        % Initial height based on number of folders
        nFolders = length(app.drawer_list);
        initialHeight = 200 + max(nFolders, 1) * 14;

        % Position next to main app if possible
        if isvalid(app.EyeFlowUIFigure)
            mainPos = app.EyeFlowUIFigure.Position;
            xPos = mainPos(1) + mainPos(3) + 20;
            yPos = mainPos(2);
        else
            xPos = 300;
            yPos = 300;
        end

        obj.Figure = uifigure('Position', [xPos, yPos, 750, initialHeight], ...
            'Color', [0.2, 0.2, 0.2], ...
            'Name', 'Folder Management - EyeFlow', ...
            'Resize', 'on', ...
            'WindowStyle', 'normal', ...
            'CloseRequestFcn', @(~, ~) obj.closeFigure());

        % Main layout: text area (top) + button panel (bottom)
        mainGrid = uigridlayout(obj.Figure, [2, 1], ...
            'ColumnWidth', {'1x'}, ...
            'RowHeight', {'1x', 'fit'}, ...
            'Padding', [10, 10, 10, 10], ...
            'BackgroundColor', [0.2, 0.2, 0.2], ...
            'RowSpacing', 10);

        % Text area for folder list
        if isempty(app.drawer_list)
            displayValue = {''};
        else
            displayValue = app.drawer_list;
        end

        obj.TextArea = uitextarea('Parent', mainGrid, ...
            'BackgroundColor', [0.2, 0.2, 0.2], ...
            'FontColor', [0.8, 0.8, 0.8], ...
            'Value', displayValue, ...
            'Editable', 'off');
        obj.TextArea.Layout.Row = 1;
        obj.TextArea.Layout.Column = 1;

        % Button panel grid (3 columns)
        btnGrid = uigridlayout(mainGrid, [5, 3], ...
            'ColumnWidth', {200, 200, 200}, ...
            'RowHeight', repmat({'fit'}, 1, 5), ...
            'Padding', [5, 5, 5, 5], ...
            'BackgroundColor', [0.2, 0.2, 0.2]);
        btnGrid.Layout.Row = 2;
        btnGrid.Layout.Column = 1;

        % Spinner label and spinner – adjust layout: place them in row1, col1-2
        uilabel(btnGrid, 'Text', 'Parallel workers:', ...
            'FontColor', [1 1 1], ...
            'HorizontalAlignment', 'left');
        obj.WorkerSpinner = uispinner(btnGrid, ...
            'Limits', [0, parcluster('local').NumWorkers], ...
            'Value', min(10, floor(parcluster('local').NumWorkers / 2)), ...
            'FontColor', [1 1 1], ...
            'BackgroundColor', [0.3 0.3 0.3]);
        obj.WorkerSpinner.Layout.Row = 1;
        obj.WorkerSpinner.Layout.Column = 2;

        % Common button style
        bkgColor = [0.5, 0.5, 0.5];
        fontColor = [1, 1, 1];
        renderColor = [0.2, 0.6, 0.2];

        % Column 1
        obj.makeButton(btnGrid, 2, 1, 'Select Folder', @(~, ~) obj.selectFolder(), bkgColor, fontColor);
        obj.makeButton(btnGrid, 3, 1, 'Select Entire Session', @(~, ~) obj.selectEntireFolder(), bkgColor, fontColor);
        obj.makeButton(btnGrid, 4, 1, 'Select Latest HD', @(~, ~) obj.selectLatestHDFolder(), bkgColor, fontColor);
        obj.makeButton(btnGrid, 5, 1, 'Clear List', @(~, ~) obj.clearList(), bkgColor, fontColor);

        % Column 2
        obj.makeButton(btnGrid, 2, 2, 'Load from txt', @(~, ~) obj.loadFromTxt(), bkgColor, fontColor);
        obj.makeButton(btnGrid, 3, 2, 'Save to txt', @(~, ~) obj.saveToTxt(), bkgColor, fontColor);
        obj.makeButton(btnGrid, 4, 2, 'Clear Parameters', @(~, ~) obj.clearParameters(), bkgColor, fontColor);
        obj.makeButton(btnGrid, 5, 2, 'Import Parameter', @(~, ~) obj.importParameter(), bkgColor, fontColor);

        % Column 3
        obj.makeButton(btnGrid, 2, 3, 'RENDER', @(~, ~) obj.render(), renderColor, fontColor);
        obj.makeButton(btnGrid, 3, 3, 'Show Results', @(~, ~) obj.showResults(), bkgColor, fontColor);
    end

    function btn = makeButton(~, parent, row, col, label, cb, bgColor, fontColor)
        btn = uibutton(parent, 'push', ...
            'BackgroundColor', bgColor, ...
            'FontColor', fontColor, ...
            'Text', label, ...
            'ButtonPushedFcn', cb);
        btn.Layout.Row = row;
        btn.Layout.Column = col;
    end

    function updateDisplay(obj)
        app = obj.MainApp;

        if isempty(app.drawer_list)
            obj.TextArea.Value = {''};
        else
            obj.TextArea.Value = app.drawer_list;
        end

        % Adjust figure height
        n = length(app.drawer_list);
        newHeight = 200 + max(n, 1) * 14;
        obj.Figure.Position(4) = newHeight;
    end

    function closeFigure(obj)
        delete(obj.Figure);
        delete(obj); % object destroyed
    end

    % -----------------------------------------------------------------
    % Callbacks
    % -----------------------------------------------------------------
    function selectFolder(obj)
        app = obj.MainApp;

        if ~isempty(app.drawer_list)
            startPath = app.drawer_list{end};
        else
            startPath = pwd;
        end

        selPath = uigetdir(startPath, 'Select a measurement folder');
        if selPath == 0, return; end
        app.drawer_list{end + 1} = selPath;
        obj.updateDisplay();
    end

    function selectEntireFolder(obj)
        app = obj.MainApp;

        if ~isempty(app.drawer_list)
            startPath = app.drawer_list{end};
        else
            startPath = pwd;
        end

        sessionPath = uigetdir(startPath, 'Select the session (parent) folder');
        if sessionPath == 0, return; end

        % List all subfolders
        d = dir(sessionPath);
        d = d([d.isdir]);
        d = d(~ismember({d.name}, {'.', '..'}));

        added = 0;

        for i = 1:length(d)
            folderName = d(i).name;
            % Keep only folders containing '_HD_' or '_HW_' (your naming convention)
            if contains(folderName, '_HD_') || contains(folderName, '_HW_')
                app.drawer_list{end + 1} = fullfile(sessionPath, folderName);
                added = added + 1;
            end

        end

        fprintf('Added %d measurement folders from session.\n', added);
        obj.updateDisplay();
    end

    function clearList(obj)
        obj.MainApp.drawer_list = {};
        obj.updateDisplay();
    end

    function loadFromTxt(obj)
        app = obj.MainApp;
        [file, path] = uigetfile('*.txt', 'Load folder list');
        if isequal(file, 0), return; end
        lines = readlines(fullfile(path, file));

        for i = 1:numel(lines)

            if ~isempty(lines(i)) && isfolder(lines(i))
                app.drawer_list{end + 1} = char(lines(i));
            end

        end

        obj.updateDisplay();
    end

    function saveToTxt(obj)
        app = obj.MainApp;

        if isempty(app.drawer_list)
            uialert(obj.Figure, 'List is empty.', 'Nothing to save');
            return;
        end

        [file, path] = uiputfile('*.txt', 'Save folder list as');
        if isequal(file, 0), return; end
        writelines(app.drawer_list, fullfile(path, file));
    end

    function clearParameters(obj)
        app = obj.MainApp;

        if isempty(app.drawer_list)
            uialert(obj.Figure, 'List is empty.', 'No folders');
            return;
        end

        choice = uiconfirm(obj.Figure, ...
            'Delete all JSON parameter files from the listed folders?', ...
            'Confirm Clear', 'Options', {'Yes', 'Cancel'}, ...
            'DefaultOption', 2, 'CancelOption', 2);
        if ~strcmp(choice, 'Yes'), return; end

        nDeleted = 0;

        for i = 1:length(app.drawer_list)
            jsonDir = fullfile(app.drawer_list{i}, 'eyeflow', 'json');

            if isfolder(jsonDir)
                files = dir(fullfile(jsonDir, '*.json'));

                for j = 1:length(files)
                    delete(fullfile(jsonDir, files(j).name));
                    nDeleted = nDeleted + 1;
                end

            end

        end

        fprintf('Deleted %d JSON parameter files.\n', nDeleted);
    end

    function importParameter(obj)
        app = obj.MainApp;

        if isempty(app.drawer_list)
            uialert(obj.Figure, 'List is empty.', 'No folders');
            return;
        end

        [file, path] = uigetfile('*.json', 'Select parameter JSON file to copy');
        if isequal(file, 0), return; end
        srcFile = fullfile(path, file);

        for i = 1:length(app.drawer_list)
            jsonDir = fullfile(app.drawer_list{i}, 'eyeflow', 'json');

            if ~isfolder(jsonDir)
                mkdir(jsonDir);
            end

            % Determine next available index
            existing = dir(fullfile(jsonDir, 'input_EF_params_*.json'));
            maxIdx = 0;

            for k = 1:length(existing)
                tok = regexp(existing(k).name, 'input_EF_params_(\d+)\.json', 'tokens');

                if ~isempty(tok)
                    idx = str2double(tok{1}{1});
                    if idx > maxIdx, maxIdx = idx; end
                end

            end

            destFile = fullfile(jsonDir, sprintf('input_EF_params_%d.json', maxIdx + 1));
            copyfile(srcFile, destFile);
        end

        fprintf('Parameter file imported to %d folders.\n', length(app.drawer_list));
    end

    function render(obj)
        app = obj.MainApp;
        folders = app.drawer_list;

        if isempty(folders)
            uialert(obj.Figure, 'No folders in list.', 'Empty List');
            return;
        end

        % Save current app state (checkboxes, etc.)
        savedSeg = app.segmentationCheckBox.Value;
        savedBFA = app.bloodFlowAnalysisCheckBox.Value;
        savedPWV = app.pulseVelocityCheckBox.Value;
        savedCSA = app.generateCrossSectionSignalsCheckBox.Value;
        savedExp = app.exportCrossSectionResultsCheckBox.Value;
        savedSpec = app.spectralAnalysisCheckBox.Value;

        % Set parallel workers from spinner
        nWorkers = obj.WorkerSpinner.Value;
        app.NumberofWorkersSpinner.Value = nWorkers;

        errors = cell(numel(folders), 1);
        errFolders = cell(numel(folders), 1);

        % Create a waitbar or just use command window output
        fprintf('\n========== BATCH PROCESSING STARTED ==========\n');

        for i = 1:length(folders)
            folder = folders{i};
            fprintf('\n--- Processing folder %d/%d: %s ---\n', i, length(folders), folder);
            tStart = tic;

            try
                app.Load(folder);
                err = app.ExecuteButtonPushed(); % This returns any error that occurred

                if ~isempty(err)
                    errors{i} = err;
                    errFolders{i} = folder;
                end

            catch ME
                errors{i} = ME;
                errFolders{i} = folder;
            end

            app.ClearButtonPushed();
            fprintf('  Elapsed: %.1f s\n', toc(tStart));
        end

        % Restore checkboxes
        app.segmentationCheckBox.Value = savedSeg;
        app.bloodFlowAnalysisCheckBox.Value = savedBFA;
        app.pulseVelocityCheckBox.Value = savedPWV;
        app.generateCrossSectionSignalsCheckBox.Value = savedCSA;
        app.exportCrossSectionResultsCheckBox.Value = savedExp;
        app.spectralAnalysisCheckBox.Value = savedSpec;

        % Summary
        fprintf('\n========== BATCH PROCESSING COMPLETED ==========\n');

        if isempty(errors)
            fprintf('All %d folders processed successfully.\n', length(folders));
        else
            fprintf(2, 'Errors occurred in %d folders:\n', length(errors));

            for e = 1:length(errors)
                fprintf(2, '  [%d] %s\n', e, errFolders{e});
                fprintf(2, '      %s\n', errors{e}.message);
            end

        end

    end

    function showResults(obj)
        app = obj.MainApp;

        if isempty(app.drawer_list)
            uialert(obj.Figure, 'List is empty.', 'No folders');
            return;
        end

        % Create an output directory inside the first folder (or ask user)
        outDir = fullfile(app.drawer_list{1}, 'Multiple_Results');

        if ~isfolder(outDir)
            mkdir(outDir);
        end

        fprintf('Generating multi-folder report in: %s\n', outDir);
        % Call your existing ShowOutputs function (must be on path)
        try
            ShowOutputs(app.drawer_list, outDir);
            fprintf('Report generation complete.\n');
        catch ME
            uialert(obj.Figure, sprintf('Error generating report:\n%s', ME.message), 'Report Error');
        end

    end

    function selectLatestHDFolder(obj)
        app = obj.MainApp;

        if ~isempty(app.drawer_list)
            startPath = app.drawer_list{end};
        else
            startPath = pwd;
        end

        sessionPath = uigetdir(startPath, 'Select the session (parent) folder');

        if sessionPath == 0
            return;
        end

        % Get all subfolders
        d = dir(sessionPath);
        d = d([d.isdir]);
        d = d(~ismember({d.name}, {'.', '..'}));

        % Map: baseName -> {maxNum, fullPath}
        latestMap = containers.Map('KeyType', 'char', 'ValueType', 'any');

        for i = 1:length(d)
            name = d(i).name;
            % Look for pattern: <base>_HD_<number>
            tokens = regexp(name, '^(.*)_HD_(\d+)$', 'tokens');

            if ~isempty(tokens)
                base = tokens{1}{1};
                num = str2double(tokens{1}{2});
                fullPath = fullfile(sessionPath, name);

                if ~isKey(latestMap, base) || num > latestMap(base).num
                    latestMap(base) = struct('num', num, 'path', fullPath);
                end

            end

        end

        if isempty(latestMap)
            uialert(obj.Figure, 'No folders matching *_HD_<number> found.', 'No HD Folders');
            return;
        end

        % Add each latest folder to the list
        keys = latestMap.keys;
        added = 0;

        for k = 1:length(keys)
            app.drawer_list{end + 1} = latestMap(keys{k}).path;
            added = added + 1;
            fprintf('Added latest HD folder: %s (suffix %d)\n', latestMap(keys{k}).path, latestMap(keys{k}).num);
        end

        fprintf('Added %d latest HD folders from session.\n', added);
        obj.updateDisplay();
    end

end

end
