function writeNumericToHDF5(path, datasetPath, value)
if ~isempty(value)

    % Handle for complex values (split them)
    if ~isreal(value)
        % Recursively
        writeNumericToHDF5(path, datasetPath + "_real", real(value));
        writeNumericToHDF5(path, datasetPath + "_imag", imag(value));
        return;
    end

    if islogical(value) % Hdf5 matlab does not support
        value = uint8(value);
    end

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
        value = int8(rescale(value) * 128);
    end

    h5create(path, datasetPath, size(value), 'Datatype', class(value), 'Chunksize', size(value), 'Deflate', 6);
    h5write(path, datasetPath, value);

    if ~isempty(mini)
        h5writeatt(path, datasetPath,"minimum",mini);
        h5writeatt(path, datasetPath,"maximum",MM);
    end
end
end