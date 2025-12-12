function plot_tagged_png_comparison(folder_list_source, tag)
% PLOT_TAGGED_PNG_COMPARISON Plot side-by-side PNG files containing a specific tag
%
% Inputs:
%   folder_list_source - Either a cell array of folder paths or a path to a text file
%                       containing a list of folder paths (one per line)
%   tag - String tag to search for in PNG filenames
%
% The function searches in each folder's /png subdirectory (up to 3 levels deep)
% for PNG files containing the specified tag, and plots them side by side
% with folder names as titles.

    % Read folder list from source
    if ischar(folder_list_source) || isstring(folder_list_source)
        % Input is a text file path
        folders = read_folder_list_from_file(folder_list_source);
    else
        % Input is a cell array of folder paths
        folders = folder_list_source;
    end
    
    % Validate folders exist
    valid_folders = {};
    for i = 1:length(folders)
        if exist(folders{i}, 'dir')
            valid_folders{end+1} = folders{i};
        else
            warning('Folder does not exist: %s', folders{i});
        end
    end
    
    if isempty(valid_folders)
        error('No valid folders found');
    end
    
   % --- PART 1: Find PNG files for each folder (MODIFIED) ---
    % png_files will now be a cell array of cell arrays, where each inner 
    % cell array contains all found file paths for a folder.
    png_file_groups = cell(1, length(valid_folders));
    
    for i = 1:length(valid_folders)
        folder_path = valid_folders{i};
        
        % Search for PNG files in folder/png subdirectory
        png_dir = fullfile(folder_path, 'png');
        if ~exist(png_dir, 'dir')
            warning('PNG subdirectory not found in: %s', folder_path);
            continue;
        end
        
        % Find PNG files with tag (up to 3 levels deep)
        found_files = find_png_files_with_tag(png_dir, tag, 3, 0);
        
        if isempty(found_files)
            warning('No PNG files containing tag "%s" found in: %s', tag, png_dir);
            png_file_groups{i} = {}; % Store an empty cell array
        else
            % Store ALL found files for this folder
            png_file_groups{i} = found_files; 
        end
    end
    
    % --- PART 2: Flatten the list for plotting ---
    
    % Prepare the final, flat lists for plotting
    all_png_files = {};
    all_plot_names = {};
    
    for i = 1:length(valid_folders)
        folder_path = valid_folders{i};
        found_files = png_file_groups{i};
        
        if isempty(found_files)
            continue;
        end
        
        % Extract folder name
        [~, folder_name, ~] = fileparts(folder_path);

        % Append each file to the flat list
        for j = 1:length(found_files)
            all_png_files{end+1} = found_files{j};
            % Name the plot with the folder and a sequence number if multiple exist
            if length(found_files) > 1
                all_plot_names{end+1} = sprintf('%s (%d/%d)', folder_name, j, length(found_files));
            else
                all_plot_names{end+1} = folder_name;
            end
        end
    end

    % --- PART 1 & 2: Data Finding and Flattening (RETAINED) ---

% ... [Your code for finding and flattening the file paths (all_png_files) 
%      and titles (all_plot_names) goes here] ...

    if isempty(all_png_files)
        error('No PNG files found containing tag "%s" in any folder', tag);
    end
    
% --- PART 3: Create Square-Like Figure with Tiled Layout (NEW LOGIC) ---
    
    N = length(all_png_files); % Total number of subplots needed
    
    % Calculate the optimal grid size (M rows x K columns)
    % M = rows (sqrt(N) rounded up)
    M = ceil(sqrt(N));
    % K = columns (ceil(N / M))
    K = ceil(N / M);
    
    % Define base figure size
    base_unit = 200; 
    
   % Create figure
    figure('Name', sprintf('PNG Files with Tag: "%s" (%dx%d grid)', tag, M, K), ...
           'NumberTitle', 'off', ...
           'Position', [100, 100, K * base_unit, M * base_unit + 50], ... % Adjust size
           'Visible', 'on');
           
    % **KEY STEP: Create the Tiled Layout Manager**
    % 'TileSpacing', 'none' and 'Padding', 'none' ensure maximum tightness.
    t = tiledlayout(M, K, 'TileSpacing', 'none', 'Padding', 'none');
    
    % Set an overall title for the figure if desired (optional)
    t.Title.String = sprintf('PNG Files for Tag: "%s"', tag);
    t.Title.Interpreter = 'none';
    
    % Loop through all images and plot them
    for i = 1:N
        % **KEY STEP: Get the next tile (subplot) to draw on**
        ax = nexttile;
    
        try
            img = imread(all_png_files{i});
            imshow(img, 'Parent', ax);
            
            % Set title using the generated plot name
            title(ax, all_plot_names{i}, 'Interpreter', 'none', 'FontSize', 8); 
        catch ME
            warning('Error reading image: %s. Error: %s', all_png_files{i}, ME.message);
            
            % Create empty plot with error message on failure
            imshow(zeros(100, 100, 3), 'Parent', ax);
            title(ax, sprintf('%s\n(Error)', all_plot_names{i}), 'Interpreter', 'none', 'Color', 'red', 'FontSize', 8);
        end
    
        % --- TIGHT PLOT ADJUSTMENTS ---
        
        % Remove ticks and coloring for a cleaner display
        ax.XTick = [];
        ax.YTick = [];
        ax.XColor = 'none';
        ax.YColor = 'none';
        
        % In tiledlayout, the 'TileSpacing' and 'Padding' settings 
        % already handle the majority of the tightness.
    end
    
    % Add overall title
    sgtitle(sprintf('PNG Files Containing Tag: "%s"', tag), ...
        'FontSize', 12, 'FontWeight', 'bold');

    
    fprintf('Successfully plotted %d images with tag "%s"\n', length(all_png_files), tag);
end

function folders = read_folder_list_from_file(txt_file_path)
% Read folder list from text file
    try
        fid = fopen(txt_file_path, 'r');
        if fid == -1
            error('Could not open file: %s', txt_file_path);
        end
        
        folders = {};
        while ~feof(fid)
            line = strtrim(fgetl(fid));
            if ~isempty(line) && line(1) ~= '#' % Skip empty lines and comments
                folders{end+1} = line;
            end
        end
        
        fclose(fid);
    catch ME
        if fid ~= -1
            fclose(fid);
        end
        rethrow(ME);
    end
end

function png_files = find_png_files_with_tag(root_dir, tag, max_depth, current_depth)
% Recursively search for PNG files containing the tag (up to max_depth levels)
    png_files = {};
    
    if current_depth > max_depth
        return;
    end
    
    % Search for PNG files in current directory
    pattern = fullfile(root_dir, sprintf('*%s*.png', tag));
    files = dir(pattern);
    
    for i = 1:length(files)
        if ~files(i).isdir
            png_files{end+1} = fullfile(root_dir, files(i).name);
        end
    end
    
    % Recursively search subdirectories if we haven't reached max depth
    if current_depth < max_depth
        subdirs = dir(root_dir);
        for i = 1:length(subdirs)
            if subdirs(i).isdir && ~strcmp(subdirs(i).name, '.') && ~strcmp(subdirs(i).name, '..')
                subdir_path = fullfile(root_dir, subdirs(i).name);
                png_files = [png_files, find_png_files_with_tag(subdir_path, tag, max_depth, current_depth + 1)];
            end
        end
    end
end