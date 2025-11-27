function writeNumericToHDF5(path, datasetPath, value)
if ~isempty(value)

    if ~isreal(value)
        warning("Complexe values should be handled before call");
        return;
    end

    if islogical(value) % Hdf5 matlab does not support
        value = uint8(value);
    end

    % Determine byte size of the array

    % typeSizes = struct('double',8,'single',4,'int64',8,'uint64',8,'int32',4,'uint32',4,...
    %                    'int16',2,'uint16',2,'int8',1,'uint8',1,'logical',1);
    % typeName = class(value);
    % bytesize = numel(value) * typeSizes.(typeName);
    stored_value = whos('value');
    bytesize = stored_value.bytes;
    original_class = stored_value.class;

    % Threshold for precision reduction
    threshold = 1e6; % 1 MBytes

    % Reduce precision if numeric type allows it
    isReduced = false;
    if bytesize > threshold
        isReduced = true;
        mini = min(value(:));
        MM = max(value(:));
        value = uint8(rescale(value) * 255);
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