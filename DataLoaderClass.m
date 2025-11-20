classdef DataLoaderClass < handle
% Handles all data loading operations

properties
    directory char
    filenames char
    M0 single
    M1 single
    M2 single
    SH single
end

methods

    function obj = DataLoaderClass(path)
        obj.directory = obj.validatePath(path);
        obj.filenames = obj.extractFilename(path);

        if ~isfolder(obj.directory)
            obj.loadHoloFile();
        else
            obj.loadFromRawDirectory();
        end

        if ~any(obj.M0)
            error('Data loading failed. Please check the input file.');
        end

    end

    function directory = validatePath(~, path)

        if ~isfolder(path)
            [filepath, name, ~] = fileparts(path);
            resultDir = fullfile(filepath, name);

            if ~isfolder(resultDir)
                mkdir(resultDir);
            end

            directory = resultDir;
        else
            directory = path;

            if ~isfolder(fullfile(path, "raw"))
                error('No raw file found at: %s\nPlease check folder path and filename.', path);
            end

        end

    end

    function filename = extractFilename(~, path)
        if ~isfolder(path)
            [~, name, ~] = fileparts(path);
            filename = name;
        else
            tmp_idx = regexp(path, '\');
            if iscell(tmp_idx), tmp_idx = tmp_idx{1}; end
            filename = extractBetween(path, tmp_idx(end-1)+1, tmp_idx(end)-1);
        end
    end

    function loadHoloFile(obj)
        fprintf('Reading moments in: %s.holo\n', obj.directory);
        [videoM0, videoM1, videoM2] = readMoments(strcat(obj.directory, '.holo'));
        readMomentsFooter(obj.directory);
        obj.M0 = pagetranspose(videoM0);
        obj.M1 = pagetranspose(videoM1 / 1e3);
        obj.M2 = pagetranspose(videoM2 / 1e6);
    end

    function loadFromRawDirectory(obj)
        dir_path_raw = fullfile(obj.directory, 'raw');
        h5_files = dir(fullfile(dir_path_raw, '*.h5'));
        raw_files = dir(fullfile(dir_path_raw, '*.raw'));

        if ~isempty(h5_files)
            obj.filereadHDF5();
        elseif ~isempty(raw_files)
            obj.filereadRaw();
        else
            error('No data file was found in the folder: %s', dir_path_raw);
        end

    end

    function delete(obj)
        obj.M0 = [];
        obj.M1 = [];
        obj.M2 = [];
        obj.SH = [];
    end

end

methods (Access = private)

    function filereadRaw(obj)
        readRaw(obj);

    end

    function filereadHDF5(obj)
        readHDF5(obj);
    end

end

end
