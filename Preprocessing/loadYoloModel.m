function [yoloModel] = loadYoloModel(model_name, hugging_face_repo)
%   loadEyeDiaphragmSegmentation - Loads YOLO model for eye diaphragm segmentation
%   Inputs:
%       model_name - The name of the model (how it is stored)
%       hugging_face_repo - The repo name where the model is stored on Hugging Face
%   Outputs:
%       yoloModel - Loaded YOLO model

model_path = sprintf('Models/%s.pt', model_name);

if ~isfile(model_path)
    % Download the model from Hugging Face
    url = sprintf('https://huggingface.co/DigitalHolography/%s/resolve/main/%s.pt?download=true', hugging_face_repo, model_name);
    fprintf("Downloading %s YOLO model in %s at %s ...\n", model_name, model_path, url);
    tic;
    websave(model_path, url);
    fprintf("Finished downloading %s in %s. It took: %.2fs\n", model_name, model_path, toc);
end

% Import the YOLO lib and the python environment
importYoloLib();

% Import the YOLO model using a custom Python function
yoloModel = py.yolo.import_model(model_path);

end
