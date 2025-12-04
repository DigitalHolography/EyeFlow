% install_toolboxes.m
% Script to install key MathWorks toolboxes and the ONNX converter add-on
%
% Run this in MATLAB (R2020b or later) with an active internet connection
% and appropriate license access.

toolboxes = {
    'Curve Fitting Toolbox'
    'Deep Learning Toolbox'
    'Image Processing Toolbox'
    'MATLAB'
    'Optimization Toolbox'
    'Parallel Computing Toolbox'
    'Signal Processing Toolbox'
    'Statistics and Machine Learning Toolbox'
    'Wavelet Toolbox'
};

fprintf('=== Checking for installed toolboxes ===\n');
installed = ver;
installedNames = {installed.Name};

for k = 1:numel(toolboxes)
    name = toolboxes{k};
    if ismember(name, installedNames)
        fprintf('[✔] %s is already installed.\n', name);
    else
        fprintf('[⏳] Installing %s...\n', name);
        try
            matlab.addons.toolbox.installToolbox(name);
            fprintf('[+] Installed %s successfully.\n', name);
        catch ME
            fprintf('[⚠] Could not install %s: %s\n', name, ME.message);
        end
    end
end

fprintf('\n=== Installing Deep Learning Toolbox Converter for ONNX Model Format ===\n');
try
    % Install the ONNX converter add-on directly from Add-On Explorer
    addonName = 'Deep Learning Toolbox Converter for ONNX Model Format';
    matlab.addons.toolbox.installToolbox(addonName);
    fprintf('[+] Installed %s successfully.\n', addonName);
catch
    fprintf('[⏳] Attempting to install from Add-On Explorer...\n');
    try
        % Search and install by name from the MathWorks Add-On Explorer
        matlab.addons.installAddonFromMAL(addonName);
        fprintf('[+] Installed %s successfully from Add-On Explorer.\n', addonName);
    catch ME
        fprintf('[❌] Could not install %s automatically.\nError: %s\n', addonName, ME.message);
        fprintf('Please open MATLAB Add-On Explorer manually and install "%s".\n', addonName);
    end
end

fprintf('\n=== Done! ===\n');
