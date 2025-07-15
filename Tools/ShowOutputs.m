function ShowOutputs(paths, output_dir)
% SHOWOUTPUTS Display multiple outputs from foldermanagement drawerlist
%   Shows various analysis results (segmentation, velocity, volume, etc.) 
%   for each path and exports them as montage images to output directory
%   - Missing figures are replaced with blank placeholders to maintain order
%
% Inputs:
%   paths - cell array of paths to process
%   output_dir - directory to save output images

% Validate inputs
if nargin < 2
    error('Both paths and output_dir arguments are required');
end

if ~iscell(paths)
    paths = {paths}; % Convert single path to cell array
end

if ~exist(output_dir, 'dir')
    mkdir(output_dir); % Create output directory if it doesn't exist
end

N = length(paths);

% Define all output types and their file patterns
output_types = {
    'segmentation', '_vesselMap.png', 'png/mask';
    'bvr', '_allrad_Artery_time.png', 'png/crossSectionsAnalysis';
    'Arteries_fRMS', '_f_artery_graph.png', 'png/bloodFlowVelocity';
    'ARI_velocity', '_RI_velocityArtery.png', 'png/bloodFlowVelocity';
    'histo_art_velocity', '_histogramVelocityArtery.png', 'png/bloodFlowVelocity';
    'Stroke_total_volume', '_strokeAndTotalVolume_Artery.png', 'png/crossSectionsAnalysis';
    'Vessels_velocity', '_vessels_velocity_graph.png', 'png/bloodFlowVelocity';
    'VRI_velocity', '_RI_velocityVein.png', 'png/bloodFlowVelocity';
    'A_sections', '_A_sections.png', 'png/crossSectionsAnalysis';
    'diasys_Artery', '_diasys_Artery.png', 'png/crossSectionsAnalysis';
    'diasys_Vein', '_diasys_Vein.png', 'png/crossSectionsAnalysis';
    'ArterialWaveformAnalysis_artery', '_ArterialWaveformAnalysis_v_artery.png', 'png/bloodFlowVelocity';
    'ArteriovenousPhaseDelay', '_arterial_venous_correlation.png' , 'png/bloodFlowVelocity';
};

% Initialize all path collections
for i = 1:size(output_types, 1)
    output_paths.(output_types{i,1}) = cell(1, N);
end

% Process each path
for path_idx = 1:N
    current_path = paths{path_idx};
    
    % Extract main folder name and construct expected EF folder name
    [~, main_foldername] = fileparts(current_path);
    folder_base = [char(main_foldername) '_EF'];
    ef_path = fullfile(current_path, 'eyeflow');
    
    % Skip if EyeFlow directory doesn't exist
    if ~exist(ef_path, 'dir')
        warning('EyeFlow directory not found in: %s', current_path);
        continue;
    end
    
    % Find all EF folders and get the latest one
    ef_folders = dir(fullfile(ef_path, [folder_base '_*']));
    if isempty(ef_folders)
        warning('No EF folders found in: %s', ef_path);
        continue;
    end
    
    % Extract numeric suffixes and find the maximum
    suffixes = regexp({ef_folders.name}, ['(?<=' folder_base '_)\d+$'], 'match', 'once');
    valid_idx = ~cellfun(@isempty, suffixes);
    if ~any(valid_idx)
        warning('No valid EF folders found in: %s', ef_path);
        continue;
    end
    
    max_suffix = max(str2double(suffixes(valid_idx)));
    last_folder_name = sprintf('%s_%d', folder_base, max_suffix);
    
    % Check for each output type
    for i = 1:size(output_types, 1)
        file_name = [last_folder_name output_types{i,2}];
        full_path = fullfile(ef_path, last_folder_name, output_types{i,3}, file_name);
        
        if exist(full_path, 'file')
            output_paths.(output_types{i,1}){path_idx} = full_path;
        else
            % Create a blank white placeholder image (RGB)
            placeholder_path = [];
            output_paths.(output_types{i,1}){path_idx} = placeholder_path;
        end
    end
end

% Create montages for each output type
[l, L] = bestMontageLayout(N);
figure_counter = 320;
figs_ids = [];
for i = 1:size(output_types, 1)
    type_name = output_types{i,1};
    current_paths = output_paths.(type_name);
    
    % Skip if no valid paths (unlikely due to placeholders)
    if all(cellfun(@isempty, current_paths))
        warning('No valid files (including placeholders) for output type: %s', type_name);
        continue;
    end
    
    % Create montage (missing files are replaced by placeholders)
    fi=figure(figure_counter);
    fi.Visible='on';
    montage(current_paths, 'Size', [l L]);
    exportgraphics(gca, fullfile(output_dir, [type_name '.png']));
    figs_ids = [figs_ids figure_counter];
    figure_counter = figure_counter + 1;
end
close(figs_ids);

end
