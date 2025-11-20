function [yoloModel] = loadYoloModel(model_name, hugging_face_repo)
%   loadYoloModel - Loads YOLO model
%   Inputs:
%       model_name - The name of the model (how it is stored)
%       hugging_face_repo - The repo name where the model is stored on Hugging Face
%   Outputs:
%       yoloModel - Loaded YOLO model

model_path = sprintf('Models/%s.onnx', model_name);

if ~isfile(model_path)
    % Download the model from Hugging Face
    url = sprintf('https://huggingface.co/DigitalHolography/%s/resolve/main/%s.onnx?download=true', hugging_face_repo, model_name);
    fprintf("Downloading %s model in %s at %s ...\n", model_name, model_path, url);
    tic;
    websave(model_path, url);
    fprintf("Finished downloading %s in %s. It took: %.2fs\n", model_name, model_path, toc);
end

% Import the ONNX network
try
    % Try the newer function first
    yoloModel = importNetworkFromONNX(model_path);
catch ME
    warning("Error in yolo model import : ")
    MEdisp(ME, ToolBox.EF_path);
    
    % Fall back to the older function
    yoloModel = importONNXNetwork(model_path);
end

end
