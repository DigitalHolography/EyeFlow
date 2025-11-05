function [EyeDiaphragmSegmentation] = loadEyeDiaphragmSegmentation()
%   loadEyeDiaphragmSegmentation - Loads AI networks for eye diaphragm segmentation
%   Outputs:
%       VesselSegmentationNet - Loaded vessel segmentation network

model_name = "eye_diaphragm";
model_path = 'Models\eye_diaphragm.pt';

if ~isfile(model_path)
    % Download the model from Hugging Face
    url = sprintf('https://huggingface.co/DigitalHolography/YOLO_EyeDiaphragm/blob/main/%s.pt', model_name);
    fprintf(url)
    websave(model_path, url);
end

pyenv; % Ensure Python environment is initialized

% Import the YOLO model using a custom Python function
EyeDiaphragmSegmentation = py.import_yolo_model(model_path);

end
