function writeStringToHDF5(path, datasetPath, str)
    % bytes = unicode2native(str, 'UTF-8');
    % if isempty(bytes)
    %     bytes = uint8(0);
    % end
    % h5create(path, datasetPath, size(bytes), 'Datatype', 'uint8', 'ChunkSize', size(bytes));
    h5create(path, datasetPath, [1 1], 'Datatype', 'string');
    h5write(path, datasetPath, str);
end