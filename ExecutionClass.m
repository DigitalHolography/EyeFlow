classdef ExecutionClass < handle

properties

    % Input data
    M0_ff_video double % M0 ff raw
    M0_data_video double % M0 raw
    M1_data_video double % M1 raw
    M2_data_video double % M2 raw
    SH_data_hypervideo double % SH raw

    % Preprocessed data
    f_RMS_video double % RMS sqrt(M2/M0) normalized input in kHz
    f_AVG_video double % AVG M1/M0 normalized input in kHz

    % Processing Flags
    is_preprocessed = false; % tells if the data has been preprocessed
    is_segmented = false;
    is_pulseAnalyzed = false;
    is_crossSectionAnalyzed = false;
    is_AllAnalyzed = false;

    % Checkbox flags
    flag_segmentation logical
    flag_bloodFlowVelocity_analysis logical
    flag_pulseWaveVelocity logical
    flag_crossSection_analysis logical
    flag_crossSection_figures logical
    flag_spectral_analysis logical
    flag_overwrite logical

    ToolBoxMaster ToolBoxClass

    Cache Cache
    Output Output

    displacementField % Displacement Field calculated with demons non rigid registration
    % frame by frame compared to the averaged image

    vRMS % Video of velocity map estimate in retinal vessels
    Q_results_A % Contain results from radial cross section analysis of retinal vessels
    Q_results_V
    v_video_RGB % Visual output for the velocity map estimate (arteries/red veins/blue)
    v_mean_RGB

    directory char % directory of input data (from HoloDoppler or HoloVibes)
    params_names cell % filenames of all the current input parameters ('input_EF_params.json' for example by default)
    param_name char % current filename
    filenames char % name id used for storing the measured rendered data

end

methods

    function obj = ExecutionClass(path)
        % Constructor for ExecutionClass.
        % Input: path - directory or .holo file path.

        if ~isfolder(path) % If the input path is a .holo file
            [filepath, name, ~] = fileparts(path);

            if ~isfolder(fullfile(filepath, name)) % Create a result folder
                mkdir(fullfile(filepath, name));
            end

            obj.directory = fullfile(filepath, name);
            obj.filenames = name;
        else % If input is a directory
            obj.directory = path;
            tmp_idx = regexp(path, '\');
            obj.filenames = obj.directory(tmp_idx(end - 1) + 1:end - 1);

            if ~isfolder(fullfile(path, "raw"))

                error('No raw file found at: %s\nPlease check folder path and filename.', path);

            end

        end

        % Load parameters
        obj.params_names = checkEyeFlowParamsFromJson(obj.directory);
        obj.param_name = obj.params_names{1}; % Default behavior

        % Load video data
        if ~isfolder(path) % If .holo file
            fprintf('Reading moments in: %s.holo\n', obj.directory);
            [videoM0, videoM1, videoM2] = readMoments(strcat(obj.directory, '.holo'));
            readMomentsFooter(obj.directory);
            obj.M0_ff_video = pagetranspose(improve_video(ff_correction(videoM0, 35), 0.0005, 2, 0));
            obj.M0_data_video = pagetranspose(videoM0);
            obj.M1_data_video = pagetranspose(videoM1 / 1e3); % Rescale M1
            obj.M2_data_video = pagetranspose(videoM2 / 1e6); % Rescale M2
        else
            dir_path_raw = fullfile(obj.directory, 'raw');

            % Search for all .h5 files in the folder
            h5_files = dir(fullfile(dir_path_raw, '*.h5'));
            raw_files = dir(fullfile(dir_path_raw, '*.raw'));

            if ~isempty(h5_files)
                readHDF5(obj);
            elseif ~isempty(raw_files)
                readRaw(obj);
            else
                error('No data file was found in the folder: %s', dir_path_raw);
            end

        end

        if ~any(obj.M0_data_video)
            error('Data loading failed. Please check the input file.');
        end

        obj.is_preprocessed = false;

        obj.Output = Output();
        obj.Output.initOutput();

        obj.Cache = Cache();
    end

    function preprocessData(obj)
        % Preprocess video data.

        fprintf("\n----------------------------------\n" + ...
            "Video PreProcessing\n" + ...
        "----------------------------------\n");

        % Initialize ToolBox and parameters
        ToolBox = obj.ToolBoxMaster;
        ToolBox.Output = obj.Output;
        ToolBox.Cache = obj.Cache;

        PreProcessTimer = tic;

        if any(isnan(obj.M0_data_video), 'all')
            error('NaN values found in M0 data. Please check the input file.');
        end

        % Preprocess the video data
        VideoRegistering(obj);
        VideoCropping(obj);
        VideoNormalizingLocally(obj);
        VideoResizing(obj);
        VideoNonRigidRegistering(obj);
        VideoInterpolating(obj);
        VideoRemoveOutliers(obj);

        obj.is_preprocessed = true;

        fprintf("\n----------------------------------\n" + ...
            "Preprocessing Complete\n" + ...
        "----------------------------------\n");

        fprintf("- Preprocess took : %ds\n", round(toc(PreProcessTimer)))

    end

    function analyzeData(obj, app)
        % Main routine for EyeFlow analysis.

        % Initialize ToolBox and parameters
        ToolBox = obj.ToolBoxMaster;
        params = ToolBox.getParams;
        ToolBox.Output = obj.Output;
        % ToolBox.Ref = obj; % handle to the Execution Class obj
        ToolBox.Cache = obj.Cache;

        if params.json.DebugMode
            profile off
            profile on
        end

        veins_analysis = params.veins_analysis;
        totalTime = tic;
        ToolBox.Output.add('NumFrames', size(obj.M0_data_video, 3), '', 0);
        ToolBox.Output.add('FrameRate', ToolBox.fs * 1000 / ToolBox.stride, 'Hz', 0);
        ToolBox.Output.add('InterFramePeriod', ToolBox.stride / ToolBox.fs / 1000, 's', 0);

        if ~isempty(ToolBox.record_time_stamps_us)
            tmp = ToolBox.record_time_stamps_us;
            ToolBox.Output.add('UnixTimestampFirst', tmp.first, 'µs');
            ToolBox.Output.add('UnixTimestampLast', tmp.last, 'µs');
        end

        if ~isfile(fullfile(ToolBox.path_gif, sprintf("%s_M0.gif", ToolBox.folder_name)))
            writeGifOnDisc(imresize(rescale(obj.M0_ff_video), 0.5), "M0")
        end

        % Mask Creation
        if obj.flag_segmentation
            fprintf("\n----------------------------------\n" + ...
                "Mask Creation\n" + ...
            "----------------------------------\n");
            createMasksTimer = tic;

            if ~isfolder(fullfile(ToolBox.path_png, 'mask'))
                mkdir(ToolBox.path_png, 'mask')
                mkdir(ToolBox.path_eps, 'mask')
                mkdir(fullfile(ToolBox.path_png, 'mask'), 'steps')
                mkdir(fullfile(ToolBox.path_eps, 'mask'), 'steps')
            end

            ToolBox.Cache.xy_barycenter = getBarycenter(obj.f_AVG_video);

            M0_ff_img = rescale(mean(obj.M0_ff_video, 3));

            try
                [~, diameter_x, diameter_y] = findPapilla(M0_ff_img);
            catch ME
                warning("Error while finding papilla : ")
                MEdisp(ME, ToolBox.EF_path)
                diameter_x = NaN;
                diameter_y = NaN;
            end

            createMasks(obj.M0_ff_video);
            ToolBox.Cache.papillaDiameter = mean([diameter_x, diameter_y]);

            % Visualize the segmentation result
            M0_RGB = ToolBox.Cache.M0_RGB;

            % Display the mask on the app if available
            if ~isempty(app)
                app.ImageDisplay.ImageSource = mat2gray(M0_RGB); % Rescale the image for display
                ax = ancestor(app.ImageDisplay, 'axes');
                axis(ax, 'equal');
            end

            obj.is_segmented = true;

            fprintf("- Mask Creation took: %ds\n", round(toc(createMasksTimer)));
        end

        % Pulse Analysis
        if obj.flag_bloodFlowVelocity_analysis

            fprintf("\n----------------------------------\n" + ...
                "Blood Flow Velocity Analysis\n" + ...
            "----------------------------------\n");
            pulseAnalysisTimer = tic;

            f_AVG_mean = squeeze(mean(obj.f_AVG_video, 3));
            [obj.vRMS, obj.v_video_RGB, obj.v_mean_RGB] = pulseAnalysis(obj.f_RMS_video, obj.M0_ff_video);

            if params.json.PulseAnalysis.ExtendedFlag
                extendedPulseAnalysis(obj.M0_ff_video, obj.f_RMS_video, f_AVG_mean, obj.vRMS);
            end

            axialAnalysis(obj.f_AVG_video);

            obj.is_pulseAnalyzed = true;

            fprintf("- Blood Flow Velocity Analysis took: %ds\n", round(toc(pulseAnalysisTimer)));
        end

        % Pulse Velocity Calculation
        if obj.flag_pulseWaveVelocity
            fprintf("\n----------------------------------\n" + ...
                "Pulse Velocity Calculation\n" + ...
            "----------------------------------\n");
            pulseVelocityTimer = tic;

            maskArtery = ToolBox.Cache.maskArtery;
            pulseVelocity(obj.M0_ff_video, obj.displacementField, maskArtery, 'artery');

            if veins_analysis
                maskVein = ToolBox.Cache.maskVein;
                pulseVelocity(obj.M0_ff_video, obj.displacementField, maskVein, 'vein');
            end

            time_pulsevelocity = toc(pulseVelocityTimer);
            fprintf("- Pulse Velocity Calculations took : %ds\n", round(time_pulsevelocity))
        end

        % Cross-Section Analysis
        if obj.flag_crossSection_analysis
            fprintf("\n----------------------------------\n" + ...
                "Cross-Section Analysis\n" + ...
            "----------------------------------\n");
            crossSectionAnalysisTimer = tic;

            maskArtery = ToolBox.Cache.maskArtery;
            [obj.Q_results_A] = crossSectionsAnalysis(maskArtery, 'artery', obj.vRMS, obj.M0_ff_video);

            if veins_analysis
                maskVein = ToolBox.Cache.maskVein;
                [obj.Q_results_V] = crossSectionsAnalysis(maskVein, 'vein', obj.vRMS, obj.M0_ff_video);
            end

            obj.is_crossSectionAnalyzed = true;

            fprintf("- Cross-Section Analysis took: %ds\n", round(toc(crossSectionAnalysisTimer)));
        end

        % Cross-Section Figures
        if obj.flag_crossSection_figures
            fprintf("\n----------------------------------\n" + ...
                "Cross-Section Figures\n" + ...
            "----------------------------------\n");
            crossSectionFiguresTimer = tic;

            crossSectionsFigures(obj.Q_results_A, 'artery', obj.M0_ff_video, obj.v_video_RGB, obj.v_mean_RGB);

            if veins_analysis
                crossSectionsFigures(obj.Q_results_V, 'vein', obj.M0_ff_video, obj.v_video_RGB, obj.v_mean_RGB);
                maskVessel = ToolBox.Cache.maskArtery | ToolBox.Cache.maskVein;
                sectionImageAdvanced(rescale(mean(obj.M0_ff_video, 3)), obj.Q_results_A.maskLabel, obj.Q_results_V.maskLabel, obj.Q_results_A.rejected_mask, obj.Q_results_V.rejected_mask, maskVessel);
            else
                maskArtery = ToolBox.Cache.maskArtery;
                sectionImageAdvanced(rescale(mean(obj.M0_ff_video, 3)), obj.Q_results_A.maskLabel, [], obj.Q_results_A.rejected_mask, [], maskArtery);
            end

            try

                if veins_analysis
                    combinedCrossSectionAnalysis(obj.Q_results_A, obj.Q_results_V, obj.M0_ff_video)
                end

            catch ME
                MEdisp(ME, ToolBox.EF_path)
            end

            obj.is_AllAnalyzed = true;

            fprintf("- Cross-Section Figures took: %ds\n", round(toc(crossSectionFiguresTimer)));
        end

        % Spectral Analysis
        if obj.flag_spectral_analysis && isfile(fullfile(ToolBox.EF_path, 'raw', [strcat(ToolBox.folder_name, '_SH'), '.raw']))

            fprintf("\n----------------------------------\n" + ...
                "Spectral Analysis\n" + ...
            "----------------------------------\n");
            timeSpectralAnalysis = tic;

            if ~isfolder(fullfile(ToolBox.path_png, 'spectralAnalysis'))
                mkdir(fullfile(ToolBox.path_png), 'spectralAnalysis');
                mkdir(fullfile(ToolBox.path_eps), 'spectralAnalysis');
            end

            % Spectrum Analysis
            fprintf("\n----------------------------------\n" + ...
                "Spectrum Analysis\n" + ...
            "----------------------------------\n");
            spectrumAnalysisTimer = tic;

            spectrum_analysis(obj.SH_data_hypervideo, obj.M0_ff_video);

            fprintf("- Spectrum Analysis took : %ds\n", round(toc(spectrumAnalysisTimer)))

            % Spectrogram
            fprintf("\n----------------------------------\n" + ...
                "Spectrogram\n" + ...
            "----------------------------------\n");
            spectrogramTimer = tic;

            % spectrogram(obj.maskArtery, obj.xy_barycenter, obj.SH_data_hypervideo);
            maskArtery = ToolBox.Cache.maskArtery;
            maskNeighbors = ToolBox.Cache.maskNeighbors;
            spectrum_video(obj.SH_data_hypervideo, obj.f_RMS_video, maskArtery, maskNeighbors);

            fprintf("- Spectrogram took: %ds\n", round(toc(spectrogramTimer)));

            fprintf("\n----------------------------------\n" + ...
                "Spectral Analysis Complete\n" + ...
            "----------------------------------\n");
            fprintf("Spectral Analysis timing: %ds\n", round(toc(timeSpectralAnalysis)));
        end

        if obj.is_pulseAnalyzed && veins_analysis

            try
                generateA4Report()
            catch ME
                MEdisp(ME, ToolBox.EF_path)
            end

        end

        % Main Output Saving

        ToolBox.Output.writeJson(fullfile(ToolBox.path_json, strcat(ToolBox.folder_name, 'output.json')));
        ToolBox.Output.writeHdf5(fullfile(ToolBox.path_json, strcat(ToolBox.folder_name, 'output.h5')));

        % Final Output
        tTotal = toc(totalTime);
        fprintf("\n----------------------------------\nTotal EyeFlow timing: %ds\n", round(tTotal));
        fprintf("End Computer Time: %s\n", datetime('now', 'Format', 'yyyy/MM/dd HH:mm:ss'));

        clear ToolBox;
        diary off;
        displaySuccessMsg(1);

        if params.json.DebugMode
            profile off
            profile viewer
        end

    end

end

end
