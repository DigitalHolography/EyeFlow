function writeStringToHDF5(path, datasetPath, str)
    bytes = unicode2native(str, 'UTF-8');
    if isempty(bytes)
        bytes = uint8(0);
    end
    h5create(path, datasetPath, size(bytes), 'Datatype', 'uint8', 'ChunkSize', size(bytes));
    h5write(path, datasetPath, bytes);
end