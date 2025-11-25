function [OpticDiskDetectorNet] = loadOpticDiskNetwork()
% loadOpticDiskNetwork Loads AI networks for optic disc segmentation
%   Outputs:
%       OpticDiskNet - Loaded optic disc segmentation network

model_name = "opticdisc.onnx";
model_path = 'Models\opticdisc.onnx';

if ~isfile('Models/opticdisc.onnx')
    % Download the model from Hugging Face
    url = sprintf('https://huggingface.co/DigitalHolography/EyeFlow_OpticDiscDetectorV2/resolve/main/%s', model_name);
    fprintf(url)
    websave('Models/opticdisc.onnx', url);
end

% Import the ONNX network
try
    % Try the newer function first
    OpticDiskDetectorNet = importNetworkFromONNX(model_path);
catch
    % Fall back to the older function
    warning('off')
    OpticDiskDetectorNet = importONNXNetwork(model_path);
    warning('on')
end

end
