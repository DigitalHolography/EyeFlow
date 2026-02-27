function export_h5_string(filename, name, str)

    datasetPath = "/" + name;

    % Ensure input is a string scalar
    str = string(str);

    % Check if dataset already exists
    datasetExists = false;

    if isfile(filename)
        try
            h5info(filename, datasetPath);
            datasetExists = true;
        catch
            datasetExists = false;
        end
    end

    % Create dataset if it does not exist
    if ~datasetExists
        h5create(filename, datasetPath, 1, 'Datatype', 'string');
    end

    % Write the string
    h5write(filename, datasetPath, str);

end