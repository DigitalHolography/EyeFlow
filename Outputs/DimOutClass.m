classdef DimOutClass < handle
    % Class to hold the extra data to the h5 output

    properties
        Data
    end

    % PUBLIC FUNCTIONS
    methods

        function obj = DimOutClass()
            keys = enumeration("DimEnumClass");
            numKeys = numel(keys);

            % Map the Enum directly to the dictionary
            obj.Data = dictionary(); % a dictionary of dictionaries containing "id", data fields

            for i = 1:numKeys
                obj.Data(keys(i)) = struct();
            end
        end

        function add(obj, name, data, dimDescription, unit, ste, silenceWarn)
            arguments
                obj
                name
                data
                dimDescription (1,:) string
                unit string = ""
                ste = NaN % for legacy only, not used
                silenceWarn double = 0
            end

            type = DimEnumClass.enum_from(data);

            if ~silenceWarn && ~checkDimNames(type, dimDescription)
                warning("Wrong length of dimDescription for data !\nExpected: %i\nGot: %i", type.rank, size(dimDescription, 2));
            end

            % Check && Sanitize if 1D array
            if isvector(data) && ~isscalar(data)
                data = sanitize1DArray(data);    
            end
            
            obj.Data(type).(sanitizeFieldName(name)).data = data;
            obj.Data(type).(sanitizeFieldName(name)).attributes.dimDescription = dimDescription;
            obj.Data(type).(sanitizeFieldName(name)).attributes.unit = unit;
        end

        function add_attributes(obj, name, key, val, dim)
            arguments
                obj
                name string
                key  string
                val
                dim DimEnumClass = DimEnumClass.empty
            end

            if isempty(dim)
                dirs = enumeration("DimEnumClass");
                for i = 1:numel(dirs)
                    if isfield(obj.Data(dirs{i}), name)
                        dim = dirs{i};
                        break;
                    end
                end
                warning("field \'%s\' not found", name);
            end

            obj.Data(dim).(name).attributes.(key) = val;
        end

        % FOR NOW DECREPATED
        function add_typed(obj, name, data, type, unit, ste)
            arguments
                obj
                name
                data
                type DimEnumClass
                unit string = ""
                ste = NaN
            end

            if type.validate(data)
                obj.Data(type).(sanitizeFieldName(name)) = data;
            else
                warning("DimOutClass:add", "DimEnumType (%s) is invalid with the data: %s\nSetting Other as default !", string(type), string(name));
                obj.Data(DimEnumClass.Other).(sanitizeFieldName(name)) = data;
            end

            obj.Data(type).(sanitizeFieldName(name)) = data;

            if isnumeric(ste)
                obj.Ste.(strcat(sanitizeFieldName(name))) = ste;
            else
                obj.Ste.(strcat(sanitizeFieldName(name))) = NaN;
            end

            obj.Unit.(strcat(sanitizeFieldName(name))) = unit;
        end

        %{
        D0_scalar
        D1_array
        D2_array
        D3_array
        D4_array
        strings
        %}

        function add_D0_scalar(obj, name, data, unit, ste)
            if nargin < 4
                unit = [];
            end
            if nargin < 5
                ste = [];
            end

            add_typed(obj, name, data, DimEnumClass.D0_scalar, unit, ste)
        end

        function add_D1_array(obj, name, data, unit, ste)
            if nargin < 4
                unit = [];
            end
            if nargin < 5
                ste = [];
            end

            add_typed(obj, name, data, DimEnumClass.D1_array, unit, ste)
        end

        function add_D2_array(obj, name, data, unit, ste)
            if nargin < 4
                unit = [];
            end
            if nargin < 5
                ste = [];
            end

            add_typed(obj, name, data, DimEnumClass.D2_array, unit, ste)
        end

        function add_D3_array(obj, name, data, unit, ste)
            if nargin < 4
                unit = [];
            end
            if nargin < 5
                ste = [];
            end

            add_typed(obj, name, data, DimEnumClass.D3_array, unit, ste)
        end

        function add_D4_array(obj, name, data, unit, ste)
            if nargin < 4
                unit = [];
            end
            if nargin < 5
                ste = [];
            end

            add_typed(obj, name, data, DimEnumClass.D4_array, unit, ste)
        end

        function add_strings(obj, name, data, unit, ste)
            if nargin < 4
                unit = [];
            end
            if nargin < 5
                ste = [];
            end

            add_typed(obj, name, data, DimEnumClass.Strings, unit, ste)
        end

        function writeHdf5(obj, path)
            folders = keys(obj.Data);

            for i = 1:numel(folders)
                writeHdf5_folder(obj, path, folders(i))
            end
        end
    end


    % PRIVATE FUNCTIONS
    methods (Access = private)
        function writeHdf5_folder(obj, path, folder)

            data = obj.Data(folder);

            if isempty(data)
                return
            end

            props = fieldnames(data);

            for i = 1:numel(props)
                fieldName = props{i};
                fieldValue = data.(fieldName).data;
                fieldAttributes = data.(fieldName).attributes;

                % Build dataset path
                datasetPath = "/" + string(folder) + "/" + antiSanitizeFieldName(fieldName);

                handleSavingData(obj, path, datasetPath, fieldName, fieldValue, fieldAttributes);
            end

        end

        function handleSavingData(obj, path, datasetPath, fieldName, fieldValue, fieldAttributes)
            arguments
                obj
                path
                datasetPath
                fieldName
                fieldValue
                fieldAttributes
            end

            % --- Case 1: Numeric or image data ---
            if isnumeric(fieldValue) || islogical(fieldValue)
                % % Handle for complex values (split them)
                if ~isreal(fieldValue)
                    % Recursively
                    handleSavingData(obj, path, datasetPath + "_real", fieldName, real(fieldValue), fieldAttributes);
                    handleSavingData(obj, path, datasetPath + "_imag", fieldName, imag(fieldValue), fieldAttributes);
                    return;
                end

                writeNumericToHDF5(path, datasetPath, fieldValue);

            % --- Case 2: Structs with known fields (e.g., yvalues, label) ---
            elseif isstruct(fieldValue)
                subFields = fieldnames(fieldValue);
                for j = 1:numel(subFields)
                    datasetPath = subFields{j};
                    subValue = fieldValue.(subName);
                    subPath = strcat(datasetPath, '_', subName);
                    if isnumeric(subValue)
                        writeNumericToHDF5(path, subPath, subValue);
                    elseif ischar(subValue) || isstring(subValue)
                        h5writeatt(path, datasetPath,subName, subValue);
                    end
                end

            % --- Case 3: Strings or labels ---
            elseif ischar(fieldValue) || isstring(fieldValue)
                writeStringToHDF5(path, datasetPath, string(fieldValue));

            else
                warning('Skipping unsupported field "%s" of type %s', fieldName, class(fieldValue));
            end
            writeAttributes(path, datasetPath, fieldAttributes);
        end

    end

end

% --- Helper: write numeric dataset ---

function writeAttributes(path, datasetPath, fieldAttributes)
    arguments
        path string
        datasetPath string
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

        h5writeatt(path, datasetPath, cur_key, cur_val);
    end
end

function res = sanitize1DArray(arr)
    % Sanitize the Array to be Python friendly (1, N) -> (N, 1)
    arguments
        arr 
    end
    
    res = arr(:);
end

function res_bool = checkDimNames(type, dimDescription)
    arguments
        type
        dimDescription
    end

    if type == DimEnumClass.Other
        res_bool = true;
        return;
    end

    res_bool = size(dimDescription, 2) == type.rank;
end

function name = sanitizeFieldName(str)
% Replace invalid characters for struct fields with descriptive tags
name = regexprep(str, '/', '_slash_'); % Replace / with _slash_
name = regexprep(name, '[^a-zA-Z0-9_]', '_'); % Replace any other invalid chars
% Ensure it starts with a letter
if ~isletter(name(1))
    name = ['x' name];
end
end

function original = antiSanitizeFieldName(safeName)
% Reverse the sanitization applied by sanitizeFieldName()
% Converts things like '_slash_' → '/'
% and removes any leading 'x' added for invalid first chars.

% Step 1: Replace descriptive tokens back to original characters
original = regexprep(safeName, '_slash_', '/');

% Step 2: Undo generic substitutions (if you added more rules)
% Example: if sanitizeFieldName replaced spaces or dashes, handle them here
% original = regexprep(original, '_dash_', '-');
% original = regexprep(original, '_space_', ' ');

% Step 3: If a leading 'x' was added to make a valid name, remove it
if startsWith(original, 'x') && ~startsWith(safeName, 'x_')
    original = original(2:end);
end
end
