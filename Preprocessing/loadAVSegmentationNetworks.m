function [model_struct] = loadAVSegmentationNetworks(params)
% loadAVSegmentationNetworks Loads AI networks for vesselness and artery/vein segmentation
%   Inputs:
%       params - Parameters structure containing segmentation settings
%   Outputs:
%       VesselSegmentationNet - Loaded vesselness segmentation network


use_python = false;
extension = ".onnx";
try
    % Try using Python
    pyver = pyenv;  % If this fails, Python isn't configured
    fprintf("Python detected: %s\n", pyver.Version);
    use_python = true;
    extension = ".pt";
catch
    warning("Python not detected or not configured properly. Out of date .onnx models will be used.");
end


if params.json.Mask.AVCorrelationSegmentationNet && params.json.Mask.AVDiasysSegmentationNet
    if use_python
        model_name = "nnwnet_av_corr_diasys";
    else
        model_name = "iternet5_av_corr_diasys";
    end
elseif params.json.Mask.AVDiasysSegmentationNet
    model_name = "iternet5_av_diasys";
elseif params.json.Mask.AVCorrelationSegmentationNet
    model_name = "iternet5_av_corr";
end

model_path = getLatestModel(model_name, extension);

[~, ~, ext] = fileparts(model_path);

% AVSegmentationNet = [];

model_struct = struct();
model_struct.use_python = false;
model_struct.use_onnx   = false;

if ext == ".pt"
    model_struct.use_python = true;
    model_struct.py_model = py.torch.jit.load(model_path);
elseif ext == ".onnx"
    model_struct.use_onnx = true;
    try
        % Try the newer function first
        model_struct.onnx_model = importNetworkFromONNX(model_path);
    catch
        % Fall back to the older function
        warning('off')
        model_struct.onnx_model = importONNXNetwork(model_path);
        warning('on')
    end
end

end
