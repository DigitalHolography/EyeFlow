classdef AnalyzerClass < handle
% Handles all analysis operations

properties
    vRMS
    v_video_RGB
    v_mean_RGB
    Q_results_A
    Q_results_V
end

properties (Access = private)
    ToolBoxMaster
    Output
    Cache
end

methods

    function obj = AnalyzerClass(output, cache)
        obj.Output = output;
        obj.Cache = cache;
    end

    function performSegmentation(obj, executionObj, app)
        fprintf("\n----------------------------------\nMask Creation\n----------------------------------\n");
        createMasksTimer = tic;

        ToolBox = obj.ToolBoxMaster;

        if ~isfolder(fullfile(ToolBox.path_png, 'mask'))
            mkdir(ToolBox.path_png, 'mask')
            mkdir(ToolBox.path_eps, 'mask')
            mkdir(fullfile(ToolBox.path_png, 'mask'), 'steps')
            mkdir(fullfile(ToolBox.path_eps, 'mask'), 'steps')
        end

        ToolBox.Cache.xy_barycenter = getBarycenter(executionObj.f_AVG);
        M0_ff_img = rescale(mean(executionObj.M0_ff, 3));

        try
            [~, diameter_x, diameter_y] = findPapilla(M0_ff_img);
        catch ME
            warning("Error while finding papilla : ")
            MEdisp(ME, ToolBox.EF_path)
            diameter_x = NaN;
            diameter_y = NaN;
        end

        createMasks(executionObj.M0_ff);
        ToolBox.Cache.papillaDiameter = mean([diameter_x, diameter_y]);

        % Visualize the segmentation result
        M0_RGB = ToolBox.Cache.M0_RGB;

        % Display the mask on the app if available
        if ~isempty(app)
            app.ImageDisplay.ImageSource = mat2gray(M0_RGB);
            ax = ancestor(app.ImageDisplay, 'axes');
            axis(ax, 'equal');
        end

        executionObj.is_segmented = true;
        fprintf("- Mask Creation took: %ds\n", round(toc(createMasksTimer)));
    end

    function performPulseAnalysis(obj, executionObj)
        fprintf("\n----------------------------------\nBlood Flow Velocity Analysis\n----------------------------------\n");
        pulseAnalysisTimer = tic;

        params = obj.ToolBoxMaster.getParams;
        [obj.vRMS, obj.v_video_RGB, obj.v_mean_RGB] = pulseAnalysis(executionObj.f_RMS, executionObj.M0_ff);

        if params.json.PulseAnalysis.ExtendedFlag
            f_AVG_mean = squeeze(mean(executionObj.f_AVG, 3));
            extendedPulseAnalysis(executionObj.M0_ff, executionObj.f_RMS, f_AVG_mean, obj.vRMS);
        end

        axialAnalysis(executionObj.f_AVG);

        executionObj.is_pulseAnalyzed = true;
        fprintf("- Blood Flow Velocity Analysis took: %ds\n", round(toc(pulseAnalysisTimer)));
    end

    function performPulseVelocityAnalysis(obj, executionObj)
        fprintf("\n----------------------------------\nPulse Velocity Calculation\n----------------------------------\n");
        pulseVelocityTimer = tic;

        params = obj.ToolBoxMaster.getParams;
        maskArtery = obj.Cache.maskArtery;
        pulseVelocity(executionObj.M0_ff, executionObj.displacementField, maskArtery, 'artery');

        if params.veins_analysis
            maskVein = obj.Cache.maskVein;
            pulseVelocity(executionObj.M0_ff, executionObj.displacementField, maskVein, 'vein');
        end

        time_pulsevelocity = toc(pulseVelocityTimer);
        fprintf("- Pulse Velocity Calculations took : %ds\n", round(time_pulsevelocity))
    end

    function performCrossSectionAnalysis(obj, executionObj)
        fprintf("\n----------------------------------\nCross-Section Analysis\n----------------------------------\n");
        crossSectionAnalysisTimer = tic;

        params = obj.ToolBoxMaster.getParams;
        maskArtery = obj.Cache.maskArtery;
        [obj.Q_results_A] = crossSectionsAnalysis(maskArtery, 'artery', obj.vRMS, executionObj.M0_ff);

        if params.veins_analysis
            maskVein = obj.Cache.maskVein;
            [obj.Q_results_V] = crossSectionsAnalysis(maskVein, 'vein', obj.vRMS, executionObj.M0_ff);
        end

        executionObj.is_crossSectionAnalyzed = true;
        fprintf("- Cross-Section Analysis took: %ds\n", round(toc(crossSectionAnalysisTimer)));
    end

    function generateCrossSectionFigures(obj, executionObj)
        fprintf("\n----------------------------------\nCross-Section Figures\n----------------------------------\n");
        crossSectionFiguresTimer = tic;

        params = obj.ToolBoxMaster.getParams;
        crossSectionsFigures(obj.Q_results_A, 'artery', executionObj.M0_ff, obj.v_video_RGB, obj.v_mean_RGB);

        if params.veins_analysis
            crossSectionsFigures(obj.Q_results_V, 'vein', executionObj.M0_ff, obj.v_video_RGB, obj.v_mean_RGB);
            maskVessel = obj.Cache.maskArtery | obj.Cache.maskVein;
            sectionImageAdvanced(rescale(mean(executionObj.M0_ff, 3)), obj.Q_results_A.maskLabel, obj.Q_results_V.maskLabel, obj.Q_results_A.rejected_mask, obj.Q_results_V.rejected_mask, maskVessel);
        else
            maskArtery = obj.Cache.maskArtery;
            sectionImageAdvanced(rescale(mean(executionObj.M0_ff, 3)), obj.Q_results_A.maskLabel, [], obj.Q_results_A.rejected_mask, [], maskArtery);
        end

        try

            if params.veins_analysis
                combinedCrossSectionAnalysis(obj.Q_results_A, obj.Q_results_V, executionObj.M0_ff)
            end

        catch ME
            MEdisp(ME, obj.ToolBoxMaster.EF_path)
        end

        executionObj.is_AllAnalyzed = true;
        fprintf("- Cross-Section Figures took: %ds\n", round(toc(crossSectionFiguresTimer)));
    end

    function performSpectralAnalysis(obj, executionObj)
        fprintf("\n----------------------------------\nSpectral Analysis\n----------------------------------\n");
        timeSpectralAnalysis = tic;

        ToolBox = obj.ToolBoxMaster;

        if ~isfolder(fullfile(ToolBox.path_png, 'spectralAnalysis'))
            mkdir(fullfile(ToolBox.path_png), 'spectralAnalysis');
            mkdir(fullfile(ToolBox.path_eps), 'spectralAnalysis');
        end

        % Spectrum Analysis
        fprintf("\n----------------------------------\nSpectrum Analysis\n----------------------------------\n");
        spectrumAnalysisTimer = tic;

        spectrum_analysis(executionObj.SH, executionObj.M0_ff);
        fprintf("- Spectrum Analysis took : %ds\n", round(toc(spectrumAnalysisTimer)))

        % Spectrogram
        fprintf("\n----------------------------------\nSpectrogram\n----------------------------------\n");
        spectrogramTimer = tic;

        maskArtery = obj.Cache.maskArtery;
        maskNeighbors = obj.Cache.maskNeighbors;
        spectrum_video(executionObj.SH, executionObj.f_RMS, maskArtery, maskNeighbors);

        fprintf("- Spectrogram took: %ds\n", round(toc(spectrogramTimer)));
        fprintf("\n----------------------------------\nSpectral Analysis Complete\n----------------------------------\n");
        fprintf("Spectral Analysis timing: %ds\n", round(toc(timeSpectralAnalysis)));
    end

end

end
