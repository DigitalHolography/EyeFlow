function readHDF5(obj)

dir_path_raw = fullfile(obj.directory, 'raw');

% Search for all .h5 files in the folder
h5_files = dir(fullfile(dir_path_raw, '*.h5'));

if isempty(h5_files)
    error('No HDF5 file was found in the folder: %s', dir_path_raw);
end

% Takes the first .h5 file found
RefRawFilePath = fullfile(dir_path_raw, h5_files(1).name);

fprintf('Reading the HDF5 file: %s\n', h5_files(1).name);

try
    info = h5info(RefRawFilePath, '/');
    dataset_names = {info.Datasets.Name};

    if ismember('moment0', dataset_names)
        fprintf('Reading the M0 data\n');
        obj.M0_raw_video = squeeze(h5read(RefRawFilePath, '/moment0'));
    else
        fprintf('Warning: moment0 dataset not found\n');
    end

    if ismember('moment1', dataset_names)
        fprintf('Reading the M1 data\n');
        obj.M1_raw_video = squeeze(h5read(RefRawFilePath, '/moment1'));
    else
        fprintf('Warning: moment1 dataset not found\n');
    end

    if ismember('moment2', dataset_names)
        fprintf('Reading the M2 data\n');
        obj.M2_raw_video = squeeze(h5read(RefRawFilePath, '/moment2'));
    else
        fprintf('Warning: moment2 dataset not found\n');
    end

    if ismember('SH', dataset_names)
        fprintf('Reading the SH data\n');
        obj.SH_data_hypervideo = squeeze(h5read(RefRawFilePath, '/SH'));
    else
        fprintf('Warning: SH dataset not found\n');
    end

catch ME

    disp(['ID: ' ME.identifier])
    rethrow(ME)

end

end
