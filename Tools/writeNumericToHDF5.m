function writeNumericToHDF5(path, datasetPath, value)
if ~isempty(value) && exist(path, 'file')
    hinf = h5info(path);

    if ~isempty(hinf.Datasets) && any(strcmp(hinf.Datasets.Name, datasetPath))
        warning('Dataset %s already exists, skipping', datasetPath);
    else
        % Determine byte size of the array
        typeSizes = struct('double',8,'single',4,'int64',8,'uint64',8,'int32',4,'uint32',4,...
                           'int16',2,'uint16',2,'int8',1,'uint8',1,'logical',1);
        typeName = class(value);
        bytesize = numel(value) * typeSizes.(typeName);
    
        % Threshold for precision reduction
        threshold = 1e6; % 1 MBytes
    
        % Reduce precision if numeric type allows it
        mini = [];
        if bytesize > threshold
            mini = min(value(:));
            MM = max(value(:));
            value = int8(rescale(value) * 32767);
        end

        h5create(path, datasetPath, size(value), 'Datatype', class(value));
        h5write(path, datasetPath, value);

        if ~isempty(mini)
            h5writeatt(path, datasetPath,"minimum",mini);
            h5writeatt(path, datasetPath,"maximum",MM);
        end
    end
end
end