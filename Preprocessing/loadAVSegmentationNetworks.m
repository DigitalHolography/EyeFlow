function [model_struct] = loadAVSegmentationNetworks(params)
% loadAVSegmentationNetworks Loads AI networks for vesselness and artery/vein segmentation
%   Inputs:
%       params - Parameters structure containing segmentation settings
%   Outputs:
%       AVSegmentationNet - Loaded artery/vein segmentation network

is_pytorch_available = false;
is_cuda_available = false;

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
            is_pytorch_available = true;
            % Check if CUDA is available
            if torch.cuda.is_available()

                try
                    % Try allocating a CUDA tensor
                    test = torch.rand(int32(1)).cuda();
                    fprintf("CUDA in PyTorch is working.\n");
                    is_cuda_available = true;
                catch ME
                    warning(ME.identifier, "PyTorch CUDA allocation failed. Falling back to CPU.\n%s", ME.message);
                    fprintf("PyTorch CPU will be used.\n");
                end

            else
                fprintf("CUDA not available. PyTorch CPU will be used.\n");
            end

        catch ME
            fprintf("PyTorch import failed. Falling back to ONNX.\n");
        end

    else
        fprintf("Python version %s is not supported ; it must be >= 3.10 and < 3.13 (3.12 recommanded). Using ONNX.", pyver.Version);
    end

catch ME
    fprintf("Python not detected or not configured. Using ONNX instead.\n");
end

extension = ".onnx";

% Determine the correct model name based on params
if params.json.Mask.AVCorrelationSegmentationNet && params.json.Mask.AVDiasysSegmentationNet

    if is_pytorch_available && ~isdeployed
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

if ~isdeployed
    model_path = getLatestModel(model_name, extension);
else
    currentScriptPath = fileparts(mfilename('fullpath'));
    projectRoot = fileparts(currentScriptPath);
    model_path = fullfile(projectRoot, 'Models', model_name + '.mat');
end

% AVSegmentationNet = [];

model_struct = struct();
model_struct.use_pytorch = false;
model_struct.use_onnx = false;
model_struct.use_cuda = false;

if extension == ".pt"
    model_struct.use_pytorch = true;
    model_struct.use_cuda = is_cuda_available; % Store whether we're using CUDA

    % Load the PyTorch model
    model_struct.py_model = py.torch.jit.load(model_path, pyargs('map_location', 'cpu'));

    % Move to GPU if CUDA is available and working
    if model_struct.use_cuda
        model_struct.py_model = model_struct.py_model.cuda();
    else
        model_struct.py_model = model_struct.py_model.cpu();
    end

elseif extension == ".onnx"
    model_struct.use_onnx = true;

    if isdeployed
        net_data = load(model_path);
        % The network is expected to be the first variable in the .mat file
        f = fieldnames(net_data);
        model_struct.onnx_model = net_data.(f{1});
    else

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

end
