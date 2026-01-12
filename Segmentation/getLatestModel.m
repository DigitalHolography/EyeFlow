function model_path = getLatestModel(model_name, extension)
%GET_LATEST_MODEL Ensure the latest version of a Hugging Face model is downloaded
%
% model_name : string, name of your model on Hugging Face (e.g., "iternet5_vesselness")
% Returns the local path to the ONNX model

models_dir = 'Models';

if ~isfolder(models_dir)
    mkdir(models_dir);
end

model_path = fullfile(models_dir, model_name + extension);
version_path = fullfile(models_dir, model_name + ".version");

% URLs on Hugging Face
base_url = "https://huggingface.co/DigitalHolography/" + model_name + "/resolve/main/";
model_url = base_url + model_name;
version_url = base_url + "version.txt";

try
    latest_version = webread(version_url);
    latest_version = strtrim(latest_version); % remove whitespace/newlines
catch
    warning("Could not read version from Hugging Face. Using local file if it exists.");
    latest_version = '';
end

% Check if we need to download
need_download = ~isfile(model_path) || ~isfile(version_path);

if ~need_download
    local_version = strtrim(fileread(version_path));
    need_download = ~strcmp(latest_version, local_version);
end

if need_download
    fprintf("Downloading model %s (version: %s)...\n", model_name, latest_version);
    websave(model_path, model_url);

    % Save version locally
    fid = fopen(version_path, 'w');
    fwrite(fid, latest_version);
    fclose(fid);
else
    fprintf("Model %s is up to date (version: %s).\n", model_name, latest_version);
end

end
