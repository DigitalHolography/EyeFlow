% Define input (Models) and output directories
rootDir = fileparts(fileparts(mfilename('fullpath')));
inputDir = fullfile(rootDir, 'Models');
outputDir = inputDir;

if ~isfolder(inputDir)
    error('The folder "%s" does not exist.', inputDir);
end

% Get list of all .onnx files in the directory
onnxFiles = dir(fullfile(inputDir, '*.onnx'));

if isempty(onnxFiles)
    fprintf('No .onnx files found in "%s".\n', inputDir);
else
    fprintf('Found %d ONNX files. Starting conversion...\n', numel(onnxFiles));
end

% Loop through and convert each ONNX file
for i = 1:numel(onnxFiles)
    onnxPath = fullfile(inputDir, onnxFiles(i).name);
    [~, baseName, ~] = fileparts(onnxFiles(i).name);
    matPath = fullfile(outputDir, [baseName, '.mat']);

    if isfile(matPath)
        fprintf('MAT file for "%s" already exists. Skipping conversion.\n', onnxFiles(i).name);
        continue;
    end

    try
        fprintf('Converting "%s" ... ', onnxFiles(i).name);

        % Import ONNX model
        net = importNetworkFromONNX(onnxPath);

        % Save to .mat file
        save(matPath, 'net');
        fprintf('Done.\n');
    catch ME
        fprintf('Failed: %s\n', ME.message);
    end

end

fprintf('Conversion completed.\n');
