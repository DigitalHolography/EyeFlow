function [model_struct] = loadAVSegmentationNetworks(params)
% loadAVSegmentationNetworks Loads AI networks for vesselness and artery/vein segmentation
%   Inputs:
%       params - Parameters structure containing segmentation settings
%   Outputs:
%       AVSegmentationNet - Loaded artery/vein segmentation network

use_python = false;

try
    % Try detecting Python
    pyver = pyenv;

    fprintf("Python detected: %s\n", pyver.Version);

    % Parse version string "3.x.y"
    v = sscanf(pyver.Version, "%d.%d.%d");
    major = v(1);
    minor = v(2);

    % Check version range: 3.10 ≤ version < 3.13
    if major == 3 && minor >= 10 && minor < 13

        try
            torch = py.importlib.import_module('torch');

            % Try allocating a CUDA tensor
            test = torch.rand(int32(1)).cuda();
            fprintf("CUDA in PyTorch is working.\n");
            use_python = true;

        catch ME
            warning(ME.identifier, "PyTorch CUDA unavailable or faulty. Falling back to ONNX.\n%s", ME.message);
        end

    else
        warning("Python version %s is not supported ; it must be >= 3.10 and < 3.13 (3.12 recommanded). Using ONNX.", pyver.Version);
    end

catch ME
    warning(ME.identifier, "Python not detected or not configured. Using ONNX instead.\n%s", ME.message);
end

extension = ".onnx";

% Determine the correct model name based on params
if params.json.Mask.AVCorrelationSegmentationNet && params.json.Mask.AVDiasysSegmentationNet

    if use_python
        model_name = "nnwnet_av_corr_diasys";
        extension = ".pt";
    else
        model_name = "iternet5_av_corr_diasys";
    end

elseif params.json.Mask.AVDiasysSegmentationNet
    model_name = "iternet5_av_diasys";
elseif params.json.Mask.AVCorrelationSegmentationNet
    model_name = "iternet5_av_corr";
else
    model_struct = struct();
    return
end

model_path = getLatestModel(model_name, extension);

% AVSegmentationNet = [];

model_struct = struct();
model_struct.use_python = false;
model_struct.use_onnx = false;

if extension == ".pt"
    model_struct.use_python = true;
    model_struct.py_model = py.torch.jit.load(model_path);
elseif extension == ".onnx"
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
