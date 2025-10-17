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
    M0_ff double
    M0 double
    M1 double
    M2 double
    SH double

    % Preprocessed Data
    f_RMS double
    f_AVG double
    displacementField

    % Analysis Results
    vRMS
    v_video_RGB
    v_mean_RGB
    Q_results_A
    Q_results_V

    % Processing Flags
    is_preprocessed = false
    is_segmented = false
    is_pulseAnalyzed = false
    is_crossSectionAnalyzed = false
    is_AllAnalyzed = false

    % Checkbox flags
    flag_segmentation logical
    flag_bloodFlowVelocity_analysis logical
    flag_pulseWaveVelocity logical
    flag_crossSection_analysis logical
    flag_crossSection_figures logical
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
        obj.M0_ff = DataLoader.M0_ff;
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

        fprintf("- Video Loading took : %ds\n", round(toc(tLoading)))
    end

    function preprocessData(obj)
        % Delegate to Preprocessor
        Preprocessor = PreprocessorClass(obj.directory, obj.filenames, obj.param_name);
        Preprocessor.preprocess(obj);

        % Copy results back for backward compatibility
        obj.M0 = Preprocessor.M0;
        obj.M1 = Preprocessor.M1;
        obj.M2 = Preprocessor.M2;
        obj.M0_ff = Preprocessor.M0_ff;
        obj.f_RMS = Preprocessor.f_RMS;
        obj.f_AVG = Preprocessor.f_AVG;
        obj.displacementField = Preprocessor.displacementField;
        obj.is_preprocessed = Preprocessor.is_preprocessed;

        clear Preprocessor;
    end

    function analyzeData(obj, app)
        % Main analysis coordinator
        % Initialize output system

        ToolBox = obj.ToolBoxMaster;
        params = ToolBox.getParams;

        ToolBox.setOutput(obj.Output);
        ToolBox.setCache(obj.Cache);
        obj.Cache.createtimeVector(ToolBox, size(obj.M0, 3))

        obj.Reporter = ReporterClass(obj);

        if params.json.DebugMode
            profile off
            profile on
        end

        totalTime = tic;

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

        if obj.flag_pulseWaveVelocity
            obj.Analyzer.performPulseVelocityAnalysis(obj);
        end

        if obj.flag_crossSection_analysis
            obj.Analyzer.performCrossSectionAnalysis(obj);
            obj.Q_results_A = obj.Analyzer.Q_results_A;
            obj.Q_results_V = obj.Analyzer.Q_results_V;
        end

        if obj.flag_crossSection_figures
            obj.Analyzer.generateCrossSectionFigures(obj);
        end

        if obj.flag_spectral_analysis && ~isempty(obj.SH)
            obj.Analyzer.performSpectralAnalysis(obj);
        end

        % Generate reports and outputs
        obj.Reporter.getA4Report(obj);
        obj.Reporter.saveOutputs();
        obj.Reporter.displayFinalSummary(totalTime);

        diary off;

        if params.json.DebugMode
            profile off
            profile viewer
        end

    end

end

end
