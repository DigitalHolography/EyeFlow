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
    
    % Find PNG files for each folder
    png_files = cell(1, length(valid_folders));
    folder_names = cell(1, length(valid_folders));
    
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
            png_files{i} = '';
        elseif length(found_files) > 1
            warning('Multiple PNG files containing tag "%s" found in: %s. Using first match.', tag, png_dir);
            png_files{i} = found_files{1};
        else
            png_files{i} = found_files{1};
        end
        
        % Extract folder name for display
        [~, folder_names{i}, ~] = fileparts(folder_path);
    end
    
    % Remove folders without valid PNG files
    valid_indices = ~cellfun(@isempty, png_files);
    png_files = png_files(valid_indices);
    folder_names = folder_names(valid_indices);
    
    if isempty(png_files)
        error('No PNG files found containing tag "%s" in any folder', tag);
    end
    
   % Create figure with subplots (no margins)
    figure('Name', sprintf('PNG Files with Tag: "%s"', tag), ...
           'NumberTitle', 'off', ...
           'Position', [100, 100, 200 * length(png_files), 400], ...
           'Visible', 'on');
    
    n = length(png_files);
    
    for i = 1:n
        % Manually control subplot positions with no margins
        ax = axes('Position', [(i-1)/n, 0, 1/n, 1]); % [left bottom width height]
    
        try
            img = imread(png_files{i});
            imshow(img, 'Parent', ax);
            title(ax, folder_names{i}, 'Interpreter', 'none', 'FontSize', 10);
        catch ME
            warning('Error reading image: %s. Error: %s', png_files{i}, ME.message);
            % Create empty plot with error message
            imshow(zeros(100, 100, 3), 'Parent', ax);
            title(ax, sprintf('%s\n(Error)', folder_names{i}), 'Interpreter', 'none', 'Color', 'red');
        end
    
        % Remove any default padding/margins
        ax.Units = 'normalized';
        ax.PositionConstraint = 'outerposition';
        ax.XTick = [];
        ax.YTick = [];
        ax.XColor = 'none';
        ax.YColor = 'none';
    end
    
    % Add overall title
    sgtitle(sprintf('PNG Files Containing Tag: "%s"', tag), ...
        'FontSize', 12, 'FontWeight', 'bold');

    
    fprintf('Successfully plotted %d images with tag "%s"\n', length(png_files), tag);
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