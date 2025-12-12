function [OpticDiskDetectorNet] = loadOpticDiskNetwork()
% loadOpticDiskNetwork Loads AI networks for optic disc segmentation
%   Priority is given to loading a .mat file from the Models directory.
%   If a .mat file is not found and the application is not deployed,
%   it looks for a .onnx file or downloads it from Hugging Face.
%
%   Outputs:
%       OpticDiskDetectorNet - Loaded optic disc segmentation network

OpticDiskDetectorNet = [];
model_name = "opticdisc";
currentScriptPath = fileparts(mfilename('fullpath'));
projectRoot = fileparts(currentScriptPath);
mat_model_path = fullfile(projectRoot, 'Models', model_name + '.mat');
onnx_model_path = fullfile(projectRoot, 'Models', model_name + '.onnx');

if isfile(mat_model_path)
    fprintf('Loading .mat network: %s\n', mat_model_path);
    net_data = load(mat_model_path);
    % The network is expected to be the first variable in the .mat file
    f = fieldnames(net_data);
    OpticDiskDetectorNet = net_data.(f{1});

elseif ~isdeployed
    fprintf('No .mat network found. Looking for .onnx version.\n');

    if ~isfile(onnx_model_path)
        % Download the model from Hugging Face
        fprintf('Downloading .onnx network from Hugging Face: %s\n', model_name);
        url = sprintf('https://huggingface.co/DigitalHolography/EyeFlow_OpticDiscDetectorV2/resolve/main/%s.onnx', model_name);
        websave(onnx_model_path, url);
    end

    try
        % Try the newer import function first
        fprintf('Importing .onnx network: %s\n', onnx_model_path);
        OpticDiskDetectorNet = importNetworkFromONNX(onnx_model_path);
    catch ME
        % If it fails, display the error and try the older function
        fprintf('importNetworkFromONNX failed: %s\n', ME.message);
        fprintf('Falling back to importONNXNetwork.\n');
        warning('off')
        OpticDiskDetectorNet = importONNXNetwork(onnx_model_path);
        warning('on')
    end

end

end
