% Clean up previous build folder to avoid conflicts
if isfolder('build')

    try
        rmdir('build', 's');
        fprintf('Removed existing build directory.\n');
    catch ME
        warning('Could not remove build directory. Is a file open?');
    end

end

mkdir('build');

% Convert ONNX models to MAT before compiling
if isfile(fullfile('Tools', 'convertONNXtoMAT.m'))
    run(fullfile('Tools', 'convertONNXtoMAT.m'));
end

% Prepare Paths for Compilation
% Add all source folders to the MATLAB path so the compiler can find dependencies
addpath('BloodFlowVelocity');
addpath(fullfile('BloodFlowVelocity', 'Elastography'));
addpath('CrossSection');
addpath('Loading');
addpath('Parameters');
addpath('Preprocessing');
addpath('Scripts');
addpath('Segmentation');
addpath('SHAnalysis');
addpath('Tools');
addpath('Outputs');

% Configure Build Options
buildOpts = compiler.build.StandaloneApplicationOptions('launchBatch.m');
buildOpts.ExecutableName = 'EyeFlowBatch';
buildOpts.OutputDir = 'build';
buildOpts.Verbose = 'on';
buildOpts.TreatInputsAsNumeric = 'on';

% --- SET ICON HERE ---
% Ensure you have created 'eyeflow_logo.ico' from your png
if isfile('eyeflow_logo.png')
    buildOpts.ExecutableIcon = 'eyeflow_logo.png';
else
    warning('Icon file "eyeflow_logo.ico" not found. Executable will have default icon.');
end

% --- ADD ASSETS ---
buildOpts.AdditionalFiles = [
                             "Parameters/DefaultEyeFlowParams.json", ...
                                 "Parameters/DefaultEyeFlowParamsBatch.json", ...
                                 "eyeflow_logo.png", ...
                                 "version.txt", ...
                                 "Models/opticdisc.mat", ...
                                 "Models/iternet5_av_diasys.mat", ...
                                 "Models/iternet5_vesselness.mat"
                             ];

% Run Compilation
fprintf('Starting compilation with icon...\n');
results = compiler.build.standaloneApplication(buildOpts);
fprintf('Compilation complete.\n');
