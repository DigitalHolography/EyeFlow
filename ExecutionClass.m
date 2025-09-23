classdef ExecutionClass < handle

properties

    M0_raw_video % M0 raw
    M1_raw_video % M1 raw
    M2_raw_video % M2 raw
    M0_ff_raw_video % M0 ff raw

    M0_data_video % M0 raw
    M1_data_video % M1 raw
    M2_data_video % M2 raw
    SH_data_hypervideo % SH raw

    f_RMS_video % RMS M2/M0
    f_AVG_video % AVG M1/M0
    M0_ff_video % M0 AVI

    is_preprocessed = false; % tells if the data has been preprocessed
    is_segmented = false;
    is_pulseAnalyzed = false;
    is_crossSectionAnalyzed = false;
    is_AllAnalyzed = false;

    sysIdxList % list of frame indexes counting cardiac cycles
    diasIdx
    sysIdx % Indexes for diastole/ systole analysis
    xy_barycenter % x y position of the ONH
    papillaDiameter
    vRMS % video estimate of velocity map in retinal vessels
    Q_results_A
    Q_results_V
    v_video_RGB
    v_mean_RGB

    maskArtery
    maskVein
    maskNeighbors

    displacementField

    directory char % directory of input data (from HoloDoppler or HoloVibes)
    params_names cell % filenames of all the current input parameters ('input_EF_params.json' for example by default)
    param_name char % current filename
    filenames char % name id used for storing the measured rendered data

    flag_segmentation
    flag_bloodFlowVelocity_analysis
    flag_bloodFlowVelocity_figures
    flag_crossSection_analysis
    flag_crossSection_figures
    flag_spectral_analysis

    OverWrite logical
    ToolBoxMaster ToolBoxClass
    Cache
    Outputs
    Signals
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
            disp(['Reading moments in: ', strcat(obj.directory, '.holo')]);
            [videoM0, videoM1, videoM2] = readMoments(strcat(obj.directory, '.holo'));
            readMomentsFooter(obj.directory);
            obj.M0_ff_raw_video = pagetranspose(improve_video(ff_correction(videoM0, 35), 0.0005, 2, 0));
            obj.M0_raw_video = pagetranspose(videoM0);
            obj.M1_raw_video = pagetranspose(videoM1 / 1e3); % Rescale M1
            obj.M2_raw_video = pagetranspose(videoM2 / 1e6); % Rescale M2
        else
            dir_path_raw = fullfile(obj.directory, 'raw');

            % Search for all .h5 files in the folder
            h5_files = dir(fullfile(dir_path_raw, '*.h5'));
            raw_files = dir(fullfile(dir_path_raw, '*.raw'));

            if ~isempty(h5_files)
                obj = readHDF5(obj);
            elseif ~isempty(raw_files)
                obj = readRaw(obj);
            else
                error('No data file was found in the folder: %s', dir_path_raw);
            end

        end

        obj.is_preprocessed = false;

        obj.Outputs = Outputs();
        obj.Outputs.initOutputs();

        obj.Cache = Cache();

        obj.Signals = Signals();
        obj.Signals.initSignals();
    end

    function obj = preprocessData(obj)
        % Preprocess video data.

        fprintf("\n----------------------------------\nVideo PreProcessing\n----------------------------------\n");

        PreProcessTimer = tic;

        % Store raw video data
        obj.M0_data_video = obj.M0_raw_video;
        obj.M0_ff_video = obj.M0_ff_raw_video;
        obj.M1_data_video = obj.M1_raw_video;
        obj.M2_data_video = obj.M2_raw_video;

        % Preprocess the video data
        obj = VideoRegistering(obj);
        obj = VideoCropping(obj);
        obj = VideoNormalizingLocally(obj);
        obj = VideoResizing(obj);
        obj = VideoNonRigidRegistering(obj);
        obj = VideoInterpolating(obj);
        obj = VideoRemoveOutliers(obj);

        obj.is_preprocessed = true;
        obj.Outputs.initOutputs();

        fprintf("- Preprocess took : %ds\n", round(toc(PreProcessTimer)))

    end

    function obj = analyzeData(obj, app)
        % Main routine for EyeFlow analysis.

        % Initialize ToolBox and parameters
        ToolBox = obj.ToolBoxMaster;
        params = ToolBox.getParams;

        if params.json.DebugMode
            profile off
            profile on
        end

        veins_analysis = params.veins_analysis;
        totalTime = tic;
        ToolBox.Outputs = obj.Outputs;
        ToolBox.Signals = obj.Signals;
        ToolBox.Cache = obj.Cache;
        ToolBox.Outputs.add('NumFrames', size(obj.M0_data_video, 3), '', 0);
        ToolBox.Outputs.add('FrameRate', ToolBox.fs * 1000 / ToolBox.stride, 'Hz', 0);
        ToolBox.Outputs.add('InterFramePeriod', ToolBox.stride / ToolBox.fs / 1000, 's', 0);

        if ~isfile(fullfile(ToolBox.path_gif, sprintf("%s_M0.gif", ToolBox.folder_name)))
            writeGifOnDisc(imresize(rescale(obj.M0_ff_video), 0.5), "M0")
        end

        % Mask Creation
        if obj.flag_segmentation
            fprintf("\n----------------------------------\nMask Creation\n----------------------------------\n");
            createMasksTimer = tic;

            if ~isfolder(fullfile(ToolBox.path_png, 'mask'))
                mkdir(ToolBox.path_png, 'mask')
                mkdir(ToolBox.path_eps, 'mask')
                mkdir(fullfile(ToolBox.path_png, 'mask'), 'steps')
                mkdir(fullfile(ToolBox.path_eps, 'mask'), 'steps')
            end

            M0_ff_img = rescale(mean(obj.M0_ff_video, 3));
            cmapArtery = ToolBox.cmapArtery;
            cmapVein = ToolBox.cmapVein;
            cmapAV = ToolBox.cmapAV;

            ToolBox.Cache.list.xy_barycenter = getBarycenter(obj.f_AVG_video);
            obj.xy_barycenter = ToolBox.Cache.list.xy_barycenter;

            try
                [~, diameter_x, diameter_y] = findPapilla(M0_ff_img);
            catch E
                warning("Error while finding papilla : ")
                disp(E)
                diameter_x = NaN;
                diameter_y = NaN;
            end

            [obj.maskArtery, obj.maskVein, obj.maskNeighbors] = ...
                createMasks(obj.M0_ff_video, obj.xy_barycenter);
            obj.papillaDiameter = mean([diameter_x, diameter_y]);

            M0_Artery = setcmap(M0_ff_img, obj.maskArtery, cmapArtery);
            M0_Vein = setcmap(M0_ff_img, obj.maskVein, cmapVein);
            M0_AV = setcmap(M0_ff_img, obj.maskArtery & obj.maskVein, cmapAV);

            M0_RGB = (M0_Artery + M0_Vein) .* ~(obj.maskArtery & obj.maskVein) + M0_AV + rescale(M0_ff_img) .* ~(obj.maskArtery | obj.maskVein);

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

            fprintf("\n----------------------------------\nBlood Flow Velocity Analysis\n----------------------------------\n");
            pulseAnalysisTimer = tic;

            f_AVG_mean = squeeze(mean(obj.f_AVG_video, 3));
            [obj.vRMS, obj.sysIdxList, obj.sysIdx, obj.diasIdx, obj.v_video_RGB, obj.v_mean_RGB] = pulseAnalysis(obj.f_RMS_video, obj.M0_ff_video, ...
                obj.maskArtery, obj.maskVein, obj.maskNeighbors, obj.xy_barycenter);

            if params.json.PulseAnalysis.ExtendedFlag
                extendedPulseAnalysis(obj.M0_ff_video, obj.f_RMS_video, f_AVG_mean, obj.vRMS, obj.maskArtery, obj.maskVein, obj.xy_barycenter, obj.sysIdxList);
            end

            axialAnalysis(obj.f_AVG_video, obj.maskArtery, obj.maskVein, obj.maskNeighbors);

            obj.is_pulseAnalyzed = true;

            fprintf("- Blood Flow Velocity Analysis took: %ds\n", round(toc(pulseAnalysisTimer)));
        end

        % Pulse Velocity Analysis
        %  if obj.flag_pulseVelocity_analysis
        %     fprintf("\n----------------------------------\nPulse Velocity Calculation\n----------------------------------\n");
        %     pulseVelocityTimer = tic;

        %     pulseVelocity(obj.M0_ff_video, obj.maskArtery);

        %     time_pulsevelocity = toc(pulseVelocityTimer);
        %     fprintf("- Pulse Velocity Calculations took : %ds\n", round(time_pulsevelocity))
        % end

        % Cross-Section Analysis
        if obj.flag_crossSection_analysis
            fprintf("\n----------------------------------\nCross-Section Analysis\n----------------------------------\n");
            crossSectionAnalysisTimer = tic;

            [obj.Q_results_A] = crossSectionsAnalysis(obj.maskArtery, 'artery', obj.vRMS, obj.M0_ff_video, obj.xy_barycenter, obj.papillaDiameter, obj.sysIdx, obj.diasIdx);

            if veins_analysis
                [obj.Q_results_V] = crossSectionsAnalysis(obj.maskVein, 'vein', obj.vRMS, obj.M0_ff_video, obj.xy_barycenter, obj.papillaDiameter, obj.sysIdx, obj.diasIdx);
            end

            obj.is_crossSectionAnalyzed = true;

            fprintf("- Cross-Section Analysis took: %ds\n", round(toc(crossSectionAnalysisTimer)));
        end

        % Cross-Section Figures
        if obj.flag_crossSection_figures
            fprintf("\n----------------------------------\nCross-Section Figures\n----------------------------------\n");
            crossSectionFiguresTimer = tic;

            crossSectionsFigures(obj.Q_results_A, 'artery', obj.M0_ff_video, obj.xy_barycenter, obj.sysIdxList, obj.sysIdx, obj.diasIdx, obj.v_video_RGB, obj.v_mean_RGB);

            if veins_analysis
                crossSectionsFigures(obj.Q_results_V, 'vein', obj.M0_ff_video, obj.xy_barycenter, obj.sysIdxList, obj.sysIdx, obj.diasIdx, obj.v_video_RGB, obj.v_mean_RGB);
                sectionImageAdvanced(rescale(mean(obj.M0_ff_video, 3)), obj.Q_results_A.maskLabel, obj.Q_results_V.maskLabel, obj.Q_results_A.rejected_mask, obj.Q_results_V.rejected_mask, obj.maskArtery | obj.maskVein);
            else
                sectionImageAdvanced(rescale(mean(obj.M0_ff_video, 3)), obj.Q_results_A.maskLabel, [], obj.Q_results_A.rejected_mask, [], obj.maskArtery);
            end

            try

                if veins_analysis
                    combinedCrossSectionAnalysis(obj.Q_results_A, obj.Q_results_V, obj.M0_ff_video, obj.sysIdxList)
                end

            catch e
                disp(e)
            end

            obj.is_AllAnalyzed = true;

            fprintf("- Cross-Section Figures took: %ds\n", round(toc(crossSectionFiguresTimer)));
        end

        % Spectral Analysis
        if obj.flag_spectral_analysis && isfile(fullfile(ToolBox.EF_path, 'raw', [strcat(ToolBox.folder_name, '_SH'), '.raw']))
            fprintf("\n----------------------------------\nSpectral Analysis\n----------------------------------\n");
            timeSpectralAnalysis = tic;

            if ~isfolder(fullfile(ToolBox.path_png, 'spectralAnalysis'))
                mkdir(fullfile(ToolBox.path_png), 'spectralAnalysis');
                mkdir(fullfile(ToolBox.path_eps), 'spectralAnalysis');
            end

            % Spectrum Analysis
            fprintf("\n----------------------------------\nSpectrum Analysis\n----------------------------------\n");
            spectrumAnalysisTimer = tic;

            spectrum_analysis(obj.SH_data_hypervideo, obj.M0_ff_video);

            fprintf("- Spectrum Analysis took : %ds\n", round(toc(spectrumAnalysisTimer)))

            % Spectrogram
            fprintf("\n----------------------------------\nSpectrogram\n----------------------------------\n");
            spectrogramTimer = tic;

            % spectrogram(obj.maskArtery, obj.xy_barycenter, obj.SH_data_hypervideo);
            spectrum_video(obj.SH_data_hypervideo, obj.f_RMS_video, obj.maskArtery, obj.maskNeighbors);

            fprintf("- Spectrogram took: %ds\n", round(toc(spectrogramTimer)));
            fprintf("\n----------------------------------\nSpectral Analysis timing: %ds\n", round(toc(timeSpectralAnalysis)));
        end

        if obj.is_pulseAnalyzed && veins_analysis

            try
                generateA4Report()
            catch e
                MEdisp(e, ToolBox.EF_path)
            end

        end

        % Main Outputs Saving

        ToolBox.Outputs.writeJson(fullfile(ToolBox.path_json, strcat(ToolBox.folder_name, '_main_outputs.json')));
        ToolBox.Signals.writeJson(fullfile(ToolBox.path_json, strcat(ToolBox.folder_name, '_main_signals.json')));
        ToolBox.Outputs.writeHdf5(fullfile(ToolBox.path_json, strcat(ToolBox.folder_name, '_main_outputs.h5')));
        ToolBox.Signals.writeHdf5(fullfile(ToolBox.path_json, strcat(ToolBox.folder_name, '_main_signals.h5')));

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
