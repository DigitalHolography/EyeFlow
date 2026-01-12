classdef ExtraClass < handle
% Class to hold the extra data to the h5 output

properties
    Data
    Ste
    Unit
end

methods

    function obj = ExtraClass()
        % Constructor for the class, fills the properties with default values
        obj.Data = []; % a struct containing "id", data fields
        obj.Ste = [];
        obj.Unit = [];
    end

    function add(obj, name, data, ste, unit)
        % Method to add a new output to the class

        obj.Data.(sanitizeFieldName(name)) = data;

        if nargin > 3 && ~isempty(ste)
            obj.Ste.(strcat(sanitizeFieldName(name))) = ste;
        end

        if nargin > 4
            obj.Unit.(strcat(sanitizeFieldName(name))) = unit;
        end

    end

    function writeHdf5(obj, path)

        data = obj.Data;

        if isempty(data)
            return
        end

        props = fieldnames(data);

        for i = 1:numel(props)
            fieldName = props{i};
            fieldValue = (data.(fieldName));

            % Build dataset path
            datasetPath = strcat('/Extra/', antiSanitizeFieldName(fieldName));

            % --- Case 1: Numeric or image data ---
            if isnumeric(fieldValue) || islogical(fieldValue)
                writeNumericToHDF5(path, datasetPath, fieldValue);

                if ~isempty(obj.Ste) && isfield(obj.Ste, fieldName)
                    steValue = obj.Ste.(fieldName);
                    writeNumericToHDF5(path, strcat(datasetPath, '_ste'), steValue);
                end

                if ~isempty(obj.Unit) && isfield(obj.Unit, fieldName)
                    unitValue = obj.Unit.(fieldName);
                    h5writeatt(path, datasetPath, 'unit', unitValue);
                end

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
                        h5writeatt(path, datasetPath, subName, subValue);
                    end

                end

                % --- Case 3: Strings or labels ---
            elseif ischar(fieldValue) || isstring(fieldValue)
                writeStringToHDF5(path, datasetPath, string(fieldValue));

            else
                warning('Skipping unsupported field "%s" of type %s', fieldName, class(fieldValue));
            end

        end

    end

end

end

% --- Helper: write numeric dataset ---

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
