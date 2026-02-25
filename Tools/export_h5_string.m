function export_h5_string(filename, name, str)
dataset_path = sprintf('/%s', name);
if isfile(filename)
    info = h5info(filename);
    dataset_exists = any(strcmp({info.Datasets.Name}, name));
else
    dataset_exists = false;
end
if ~dataset_exists
    h5create(filename, dataset_path, size(str), 'Datatype', 'char');
end
h5write(filename, dataset_path, str);
end