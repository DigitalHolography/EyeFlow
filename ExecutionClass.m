classdef ExecutionClass < handle
% Main controller class for EyeFlow execution and analysis.
% Coordinates between specialized modules.

properties
    % Basic Info
    path char
    directory char
    params_names cell
    param_name char
    filenames char

    % Modules
    Analyzer AnalyzerClass
    Reporter ReporterClass

    % Input Data
    M0 single
    M1 single
    M2 single
    SH single

    % AI Networks
    AINetworks AINetworksClass

    % Preprocessed Data
    M0_ff single
    M0_ff_img single
    f_RMS single
    f_AVG single
    displacementField

    % Analysis Results
    vRMS single
    v_video_RGB uint8
    v_mean_RGB uint8
    Q_results_A
    Q_results_V
    eye_side = "none"

    % Processing Flags
    is_preprocessed = false
    is_segmented = false
    is_velocityAnalyzed = false
    is_volumeRateAnalyzed = false
    is_AllAnalyzed = false

    % Checkbox flags
    flag_segmentation logical
    flag_bloodFlowVelocity_analysis logical
    flag_pulseWaveVelocity logical
    flag_crossSection_analysis logical
    flag_crossSection_export logical
    flag_spectral_analysis logical

    % Components
    ToolBoxMaster ToolBoxClass
    Output OutputClass
    Cache CacheClass
end

methods

    function obj = ExecutionClass(path)
        % Constructor for the class, initializes all components
        tLoading = tic;
        fprintf("\n----------------------------------\n" + ...
            "Video Loading\n" + ...
        "----------------------------------\n");

        % Initialize DataLoader first
        DataLoader = DataLoaderClass(path);

        % Copy data to ExecutionClass for backward compatibility
        obj.directory = DataLoader.directory;
        obj.filenames = DataLoader.filenames;
        obj.path = path;

        % Copy data to ExecutionClass for backward compatibility
        obj.M0 = DataLoader.M0;
        obj.M1 = DataLoader.M1;
        obj.M2 = DataLoader.M2;
        obj.SH = DataLoader.SH;
        clear DataLoader;

        % Load parameters
        obj.params_names = checkEyeFlowParamsFromJson(obj.directory);
        obj.param_name = obj.params_names{1};

        % Initialize Output and Cache
        obj.Output = OutputClass();
        obj.Cache = CacheClass();

        % Initialize Modules
        obj.Analyzer = AnalyzerClass();

        Parameters_json(path, obj.param_name);

        fprintf("- Video Loading took : %ds\n", round(toc(tLoading)))
    end

    function preprocessData(obj)

        PreProcessTimer = tic;
        fprintf("\n----------------------------------\nVideo PreProcessing\n----------------------------------\n");

        % Delegate to Preprocessor
        Preprocessor = PreprocessorClass(obj.directory, obj.filenames, obj.param_name);

        % Preprocessing
        Preprocessor.preprocess(obj);

        % Copy results back for backward compatibility
        obj.M0_ff = Preprocessor.M0_ff;
        obj.M0_ff_img = squeeze(mean(obj.M0_ff, 3));
        obj.f_RMS = Preprocessor.f_RMS;
        obj.f_AVG = Preprocessor.f_AVG;
        obj.displacementField = Preprocessor.displacementField;
        obj.is_preprocessed = Preprocessor.is_preprocessed;

        % Clear Preprocessor and intermediate variables
        obj.M0 = [];
        obj.M1 = [];
        obj.M2 = [];
        clear Preprocessor;

        fprintf("- Preprocess took : %ds\n", round(toc(PreProcessTimer)))

        fprintf("\n----------------------------------\nPreprocessing Complete\n----------------------------------\n");
    end

    function analyzeData(obj, app)
        % Main analysis coordinator
        % Initialize output system

        AnalyzerTimer = tic;

        ToolBox = obj.ToolBoxMaster;
        params = ToolBox.getParams;

        ToolBox.setOutput(obj.Output);
        ToolBox.setCache(obj.Cache);

        obj.Cache.createtimeVector(ToolBox, size(obj.M0_ff, 3));

        obj.Reporter = ReporterClass(obj);

        if params.json.DebugMode
            profile off
            profile on
        end

        obj.AINetworks.updateAINetworks(params);

        % Execute eye side analysis if asked
        if params.json.Mask.EyeSideClassifierNet
            obj.eye_side = predictEyeSide(obj.AINetworks.EyeSideClassifierNet, obj.M0_ff_img);
        end

        % Execute analysis steps based on checkbox flags
        if obj.flag_segmentation
            obj.Analyzer.performSegmentation(obj, app);
        end

        if obj.flag_bloodFlowVelocity_analysis
            obj.Analyzer.performPulseAnalysis(obj);
            obj.vRMS = obj.Analyzer.vRMS;
            obj.v_video_RGB = obj.Analyzer.v_video_RGB;
            obj.v_mean_RGB = obj.Analyzer.v_mean_RGB;
        end

        try

            if obj.flag_pulseWaveVelocity
                obj.Analyzer.performPulseVelocityAnalysis(obj);
            end

        catch ME
            MEdisp(ME, '');
            warning("Pulse Wave Velocity Analysis failed");
        end

        try

            if obj.flag_crossSection_analysis
                obj.Analyzer.performCrossSectionAnalysis(obj);
                obj.Q_results_A = obj.Analyzer.Q_results_A;
                obj.Q_results_V = obj.Analyzer.Q_results_V;
            end

            if obj.flag_crossSection_export
                obj.Analyzer.generateexportCrossSectionResults(obj);
            end

        catch ME
            MEdisp(ME, "", "WARN", "Cross-Section Analysis failed");
        end

        if obj.flag_spectral_analysis && ~isempty(obj.SH)
            obj.Analyzer.performSpectralAnalysis(obj);
        end

        diary off;

        if params.json.DebugMode
            profile off
            profile viewer
        end

        fprintf("- Total Analysis took : %ds\n", round(toc(AnalyzerTimer)))

    end

    function checkData(obj)
        % Visual check of loaded/preprocessed data
        figure

        if ~isempty(obj.M0_ff)
            subplot(2, 2, 1)
            imagesc(mean(obj.M0_ff, 3))
            axis image off
            subplot(2, 2, 2)
            imagesc(mean(obj.f_AVG, 3))
            axis image off
            subplot(2, 2, 3)
            imagesc(mean(obj.f_RMS, 3))
            axis image off

            if ~isempty(obj.displacementField)
                subplot(2, 2, 4)
                imagesc(sqrt(mean(obj.displacementField .^ 2, 4)))
                axis image off
            end

        elseif ~isempty(obj.M0)
            subplot(1, 3, 1)
            imagesc(mean(obj.M0, 3))
            axis image off
            subplot(1, 3, 2)
            imagesc(mean(obj.M1, 3))
            axis image off
            subplot(1, 3, 3)
            imagesc(mean(obj.M2, 3))
            axis image off
        end

    end

end

end
