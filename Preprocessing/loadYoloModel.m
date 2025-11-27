function [yoloModel] = loadYoloModel(model_name, hugging_face_repo)
%   loadYoloModel - Loads YOLO model
%   Inputs:
%       model_name - The name of the model (how it is stored)
%       hugging_face_repo - The repo name where the model is stored on Hugging Face
%   Outputs:
%       yoloModel - Loaded YOLO model


currentScriptPath = fileparts(mfilename('fullpath'));
projectRoot = fileparts(currentScriptPath);
mat_model_path = fullfile(projectRoot, 'Models', model_name + '.mat');
onnx_model_path = fullfile(projectRoot, 'Models', model_name + '.onnx');

if isfile(mat_model_path)
    fprintf('Loading .mat network: %s\n', mat_model_path);
    net_data = load(mat_model_path);
    % The network is expected to be the first variable in the .mat file
    f = fieldnames(net_data);
    yoloModel = net_data.(f{1});
elseif ~isdeployed
    fprintf('No .mat network found. Looking for .onnx version.\n');
    if ~isfile(model_path)
        % Download the model from Hugging Face
        url = sprintf('https://huggingface.co/DigitalHolography/%s/resolve/main/%s.onnx?download=true', hugging_face_repo, model_name);
        fprintf("Downloading %s model in %s at %s ...\n", model_name, onnx_model_path, url);
        tic;
        websave(onnx_model_path, url);
        fprintf("Finished downloading %s in %s. It took: %.2fs\n", model_name, onnx_model_path, toc);
    end
    
    % Import the ONNX network
    try
        % Try the newer function first
        yoloModel = importNetworkFromONNX(onnx_model_path);
    catch ME
        warning("Error in yolo model import : ")
        MEdisp(ME, ToolBox.EF_path);
        
        % Fall back to the older function
        yoloModel = importONNXNetwork(onnx_model_path);
    end
end

end
