classdef OutputClass < handle
% Class to hold the Output of the retinal flow analysis pipeline

properties
    data
end

methods

    function obj = OutputClass()
        structTemplate = struct(...
            "value", [], ...
            "h5path", "", ...
            "metadata", [], ...
            "attributes", [] ...
            );
        obj.data = dictionary(string, structTemplate);
    end

    function add(obj, name, value, unit, ste, vars)
        arguments
            obj
            name string
            value
            unit = "" % old behavior
            ste = NaN % old behavior
            vars.h5path string = ""
            vars.unit string = ""
            vars.standard_error = NaN
            vars.keepSize = false
        end

        if nargin > 3
            vars.unit = unit;
        end

        if nargin > 4
            vars.standard_error = ste;
        end

        entry.value = value;
        entry.metadata.standard_error = vars.standard_error;
        entry.metadata.unit = vars.unit;
        entry.metadata.keepSize = vars.keepSize;
        entry.attributes = [];

        obj.data(name) = entry;

        % obj.data(name).value = value;
        % obj.data(name).metadata.standard_error = vars.standard_error;
        % obj.data(name).metadata.unit = vars.unit;

        if vars.h5path == ""
            obj.data(name).h5path = sprintf("/%s", name);
        else
            obj.data(name).h5path = (vars.h5path);
        end

    end

    function add_attribute(obj, value_name, key, val)
        arguments
            obj
            value_name string   % Name of the value
            key string          % Name of the attribute
            val
        end

        entry = obj.data(value_name);
        entry.attributes.(key) = val;
        obj.data(value_name) = entry;
    end

    function writeJson(obj, path)
        props = keys(obj.data);

        d = containers.Map('KeyType', 'char', 'ValueType', 'any');

        for i = 1:length(props)
            dic_key = props(i);
            
            if dic_key == ""
                continue;
            end

            v = obj.data(dic_key).value;

            if isscalar(v)

                if obj.data(dic_key).metadata.unit == ""
                    key = dic_key;
                else
                    key = dic_key + "_" + obj.data(dic_key).metadata.unit;
                end

                d(char(key)) = v;
            end

        end

        if ~isempty(d)
            jsonText = jsonencode(d, "PrettyPrint", true);

            fid = fopen(path, 'w');

            if fid == -1
                error('Cannot open file for writing: %s', path);
            end

            fwrite(fid, jsonText, 'char');
            fclose(fid);
        end

    end

    function writeHdf5(obj, path)

        [folder_dir, folder_name, ~] = fileparts(path);
        file_path = fullfile(folder_dir, strcat(folder_name, ".h5"));

        if isfile(file_path)
            delete(file_path)
        end

        props = keys(obj.data);

        for i = 1:length(props)
            dic_key = props(i);
            
            if dic_key == ""
                continue;
            end

            h5path = (obj.data(dic_key).h5path);
            temp = char(h5path);

            if temp(1) ~= '/'
                h5path = strcat("/", h5path);
            end

            keepSize = obj.data(dic_key).metadata.keepSize;

            if isnumeric(obj.data(dic_key).metadata.standard_error)

                if ~isnan(obj.data(dic_key).metadata.standard_error)
                    writeNumericToHDF5(file_path, strcat(h5path, "/value"), obj.data(dic_key).value, keepSize);
                    h5writeatt(file_path, strcat(h5path, "/value"), "unit", obj.data(dic_key).metadata.unit);
                    writeNumericToHDF5(file_path, strcat(h5path, "/ste"), obj.data(dic_key).metadata.standard_error, keepSize);
                    h5writeatt(file_path, strcat(h5path, "/ste"), "unit", obj.data(dic_key).metadata.unit);
                else
                    writeNumericToHDF5(file_path, strcat(h5path, "/value"), obj.data(dic_key).value, keepSize);
                    h5writeatt(file_path, strcat(h5path, "/value"), "unit", obj.data(dic_key).metadata.unit);
                end

            else
                warning_s("Non numeric ste given to save property : %s", dic_key);
                writeNumericToHDF5(file_path, h5path, obj.data(dic_key).value);
                h5writeatt(file_path, h5path, "unit", obj.data(dic_key).metadata.unit);
            end

            h5writeatt(file_path, h5path + "/value", "nameID", dic_key);

            writeAttributes(file_path, h5path, obj.data(dic_key).attributes);

        end

    end

end

end

% --- Helper: write numeric dataset ---

function writeNumericToHDF5(path, datasetPath, value, keepSize)
    arguments
        path,
        datasetPath,
        value,
        keepSize = false
    end

if ~isempty(value)

    if isnumeric(value) & ~isreal(value)
        warning("Complex values should be handled before call : %s", datasetPath);
        return;
    end

    if islogical(value) % Hdf5 matlab does not support
        value = uint8(value);
    end
    
    % Reduce precision if numeric type allows it
    isReduced = false;
    
    if ~keepSize
        stored_value = whos('value');
        bytesize = stored_value.bytes;
        original_class = stored_value.class;

        % Threshold for precision reduction
        threshold = 1e6; % 1 MBytes

        if bytesize > threshold
            isReduced = true;
            mini = min(value(:));
            MM = max(value(:));
            value = uint8(rescale(value) * 255);
        end
    end

    if isa(value, 'double')
        value = single(value);
    end

    if isfile(path)
        try
            h5info(path, datasetPath);
            datasetExists = true;
        catch
            datasetExists = false;
        end

        if datasetExists
            warning_s("Dataset '%s' already exists in '%s'. Skipping.", datasetPath, path);
            return;
        end
    end

    h5create(path, datasetPath, size(value), 'Datatype', class(value), 'Chunksize', size(value), 'Deflate', 6);
    h5write(path, datasetPath, value);

    if isReduced
        h5writeatt(path, datasetPath, "minimum", mini);
        h5writeatt(path, datasetPath, "maximum", MM);
        h5writeatt(path, datasetPath, "original_class", original_class);
    end

end

end


% ==================== Helper to write the h5 attributes ==================== %

function writeAttributes(file_path, h5path, fieldAttributes)
    arguments
        file_path string
        h5path string
        fieldAttributes struct
    end

    keys = fieldnames(fieldAttributes);

    for i = 1:numel(keys)
        cur_key = keys{i};
        cur_val = fieldAttributes.(cur_key);
        % This will avoid [strings] and "" (will replace by null)
        if isstring(cur_val) && isscalar(cur_val)
            cur_val = char(cur_val);
        end

        h5writeatt(file_path, h5path + "/value", cur_key, cur_val);
    end
end
