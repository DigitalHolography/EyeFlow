classdef AnalyzerClass < handle
% Handles all analysis operations

properties
    vRMS
    v_video_RGB
    v_mean_RGB
    Q_results_A
    Q_results_V
end

methods

    function obj = AnalyzerClass()

    end

    function performSegmentation(~, executionObj, app)
        fprintf("\n----------------------------------\n" + ...
            "Segmentation\n" + ...
        "----------------------------------\n");
        createMasksTimer = tic;

        ToolBox = getGlobalToolBox;

        if ~isfolder(fullfile(ToolBox.path_png, 'mask'))
            mkdir(ToolBox.path_png, 'mask')
        end

        if ~isfolder(fullfile(ToolBox.path_eps, 'mask'))
            mkdir(ToolBox.path_eps, 'mask')
        end

        if ~isfolder(fullfile(ToolBox.path_png, 'mask', 'steps'))
            mkdir(fullfile(ToolBox.path_png, 'mask'), 'steps')
        end

        if ~isfolder(fullfile(ToolBox.path_eps, 'mask', 'steps'))
            mkdir(fullfile(ToolBox.path_eps, 'mask'), 'steps')
        end

        ToolBox.Cache.xy_barycenter = getBarycenter(executionObj.f_AVG);
        M0_ff_img = rescale(mean(executionObj.M0_ff, 3));

        ToolBox.Output.Extra.add("M0_ff_img", M0_ff_img);

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

        params = ToolBox.getParams;
        % Artery score
        scoreA = ToolBox.Cache.scoreMaskArtery;
        if isempty(scoreA) || isnan(scoreA)
        else
            fprintf("- Mask artery quality score: %.2f\n", scoreA);
            if scoreA < 0.5
            error("AnalyzerClass:LowMaskScore", "Mask artery quality score too low (%.2f < 0.5). Aborting segmentation.", scoreA);
            end
        end
        ToolBox.Output.Extra.add("maskArtery", ToolBox.Cache.maskArtery);

        % Vein score (if veins analysis enabled)
        if params.veins_analysis
            scoreV = ToolBox.Cache.scoreMaskVein;
            if isempty(scoreV) || isnan(scoreV)
            else
                fprintf("- Mask vein quality score: %.2f\n", scoreV);
                if scoreV < 0.5
                    error("AnalyzerClass:LowMaskScore", "Mask vein quality score too low (%.2f < 0.5). Aborting segmentation.", scoreV);
                end
            end
        end
        ToolBox.Output.Extra.add("maskVein", ToolBox.Cache.maskVein);


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
        fprintf("\n----------------------------------\n" + ...
            "Blood Flow Velocity Analysis\n" + ...
        "----------------------------------\n");
        pulseAnalysisTimer = tic;

        ToolBox = getGlobalToolBox;
        params = ToolBox.getParams;
        [obj.vRMS, obj.v_video_RGB, obj.v_mean_RGB] = pulseAnalysis(executionObj.f_RMS, executionObj.M0_ff);

        if params.json.PulseAnalysis.ExtendedFlag
            f_AVG_mean = squeeze(mean(executionObj.f_AVG, 3));
            extendedPulseAnalysis(executionObj.M0_ff, executionObj.f_RMS, f_AVG_mean, obj.vRMS);
        end

        axialAnalysis(executionObj.f_AVG);

        executionObj.is_velocityAnalyzed = true;
        fprintf("- Blood Flow Velocity Analysis took: %ds\n", round(toc(pulseAnalysisTimer)));
    end

    function performPulseVelocityAnalysis(~, executionObj)
        fprintf("\n----------------------------------\n" + ...
            "Pulse Velocity Calculation\n" + ...
        "----------------------------------\n");
        pulseVelocityTimer = tic;

        ToolBox = getGlobalToolBox;
        params = ToolBox.getParams;
        maskArtery = ToolBox.Cache.maskArtery;
        pulseVelocity(executionObj.M0_ff, executionObj.displacementField, maskArtery, 'artery');

        if params.veins_analysis
            maskVein = ToolBox.Cache.maskVein;
            pulseVelocity(executionObj.M0_ff, executionObj.displacementField, maskVein, 'vein');
        end

        time_pulsevelocity = toc(pulseVelocityTimer);
        fprintf("- Pulse Velocity Calculations took : %ds\n", round(time_pulsevelocity))
    end

    function performCrossSectionAnalysis(obj, executionObj)
        fprintf("\n----------------------------------\n" + ...
            "Generate Cross-Section Signals\n" + ...
        "----------------------------------\n");
        crossSectionAnalysisTimer = tic;

        ToolBox = getGlobalToolBox;
        params = ToolBox.getParams;
        maskArtery = ToolBox.Cache.maskArtery;
        [obj.Q_results_A] = generateCrossSectionSignals(maskArtery, 'artery', obj.vRMS, executionObj.M0_ff);

        if params.veins_analysis
            maskVein = ToolBox.Cache.maskVein;
            [obj.Q_results_V] = generateCrossSectionSignals(maskVein, 'vein', obj.vRMS, executionObj.M0_ff);
        end

        executionObj.is_volumeRateAnalyzed = true;
        fprintf("- Cross-Section Signals Generation took: %ds\n", round(toc(crossSectionAnalysisTimer)));
    end

    function generateexportCrossSectionResults(obj, executionObj)
        fprintf("\n----------------------------------\n" + ...
            "Export Cross-Section Results\n" + ...
        "----------------------------------\n");
        exportCrossSectionResultsTimer = tic;

        ToolBox = getGlobalToolBox;
        params = ToolBox.getParams;
        exportCrossSectionResults(obj.Q_results_A, 'artery', executionObj.M0_ff, obj.v_video_RGB, obj.v_mean_RGB);

        if params.veins_analysis
            exportCrossSectionResults(obj.Q_results_V, 'vein', executionObj.M0_ff, obj.v_video_RGB, obj.v_mean_RGB);
            maskVessel = ToolBox.Cache.maskArtery | ToolBox.Cache.maskVein;
            sectionImageAdvanced(rescale(mean(executionObj.M0_ff, 3)), obj.Q_results_A.maskLabel, obj.Q_results_V.maskLabel, obj.Q_results_A.rejected_mask, obj.Q_results_V.rejected_mask, maskVessel);
        else
            maskArtery = ToolBox.Cache.maskArtery;
            sectionImageAdvanced(rescale(mean(executionObj.M0_ff, 3)), obj.Q_results_A.maskLabel, [], obj.Q_results_A.rejected_mask, [], maskArtery);
        end

        try

            if params.veins_analysis
                combinedCrossSectionAnalysis(obj.Q_results_A, obj.Q_results_V, executionObj.M0_ff)
            end

        catch ME
            MEdisp(ME, ToolBox.EF_path)
        end

        executionObj.is_AllAnalyzed = true;
        fprintf("- Cross-Section Results Export took: %ds\n", round(toc(exportCrossSectionResultsTimer)));
    end

    function performSpectralAnalysis(~, executionObj)
        fprintf("\n----------------------------------\n" + ...
            "Spectral Analysis\n" + ...
        "----------------------------------\n");
        timeSpectralAnalysis = tic;

        ToolBox = getGlobalToolBox;

        if ~isfolder(fullfile(ToolBox.path_png, 'spectralAnalysis'))
            mkdir(fullfile(ToolBox.path_png), 'spectralAnalysis');
            mkdir(fullfile(ToolBox.path_eps), 'spectralAnalysis');
        end

        % Spectrum Analysis
        fprintf("\n----------------------------------\n" + ...
            "Spectrum Analysis\n" + ...
        "----------------------------------\n");
        spectrumAnalysisTimer = tic;

        spectrum_analysis(executionObj.SH, executionObj.M0_ff);
        fprintf("- Spectrum Analysis took : %ds\n", round(toc(spectrumAnalysisTimer)))

        % Spectrogram
        fprintf("\n----------------------------------\n" + ...
            "Spectrogram\n" + ...
        "----------------------------------\n");
        spectrogramTimer = tic;

        maskArtery = ToolBox.Cache.maskArtery;
        maskNeighbors = ToolBox.Cache.maskNeighbors;
        spectrum_video(executionObj.SH, executionObj.f_RMS, maskArtery, maskNeighbors);

        fprintf("- Spectrogram took: %ds\n", round(toc(spectrogramTimer)));
        fprintf("\n----------------------------------\n" + ...
            "Spectral Analysis Complete\n" + ...
        "----------------------------------\n");
        fprintf("Spectral Analysis timing: %ds\n", round(toc(timeSpectralAnalysis)));
    end

end

end
