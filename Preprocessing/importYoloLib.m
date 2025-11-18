function importYoloLib()
%   Import the yolo.py python file as a library to use in MATLAB

% The folder in which yolo.py is. (by default: "Preprocessing")
lib_folder = "Preprocessing";

if ~isfolder(lib_folder)
    error('No directory called "%s"', lib_folder);
end

pyenv; % Ensure Python environment is initialized

% Add lib_folder to PATH if not already present
if count(py.sys.path, lib_folder) == 0
    py.sys.path().insert(int32(0), lib_folder);
end

mod = py.importlib.import_module('yolo');
py.importlib.reload(mod);

end