function [VesselSegmentationNet] = loadVesselSegmentationNetworks()
% loadVesselSegmentationNetworks Loads AI networks for vessel segmentation
%   Outputs:
%       VesselSegmentationNet - Loaded vessel segmentation network

model_name = "iternet5_vesselness";
model_path = 'Models\iternet5_vesselness.onnx';

if ~isfile(model_path)
    % Download the model from Hugging Face
    url = sprintf('https://huggingface.co/DigitalHolography/%s/resolve/main/%s', model_name, model_name);
    fprintf(url)
    websave(model_path, url);
end

% Import the ONNX network
try
    % Try the newer function first
    VesselSegmentationNet = importNetworkFromONNX(model_path);
catch
    % Fall back to the older function
    warning('off')
    VesselSegmentationNet = importONNXNetwork(model_path);
    warning('on')
end

fprintf("Loaded iternet5_vesselness model for vessel segmentation\n");

end
