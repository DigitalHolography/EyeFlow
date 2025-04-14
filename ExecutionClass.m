classdef ExecutionClass < handle
    
    properties
        
        M0_raw_video % M0 raw
        M1_raw_video % M1 raw
        M2_raw_video % M2 raw
        M0_ff_raw_video % M0 ff raw
        
        M0_data_video % M0 raw
        M1_data_video % M1 raw
        M2_data_video % M2 raw
        
        f_RMS_video % RMS M2/M0
        f_AVG_video % AVG M1/M0
        M0_ff_video % M0 AVI
        
        is_preprocessed = false; % tells if the data has been preprocessed
        is_segmented = false;
        is_pulseAnalyzed = false;
        is_crossSectionAnalyzed = false;
        
        sysIdxList % list of frame indexes counting cardiac cycles
        diasIdx
        sysIdx % Indexes for diastole/ systole analysis
        xy_barycenter % x y position of the ONH
        vRMS % video estimate of velocity map in retinal vessels
        Q_results_A
        Q_results_V
        
        SH_data_hypervideo % SH raw
        
        maskArtery
        maskVein
        maskNeighbors
        
        directory char % directory of input data (from HoloDoppler or HoloVibes)
        params_names cell % filenames of all the current input parameters ('InputEyeFlowParams.json' for example by default)
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
        Outputs
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
                obj = readRaw(obj);
            end
            
            obj.is_preprocessed = false;

            obj.Outputs = Outputs();
            obj.Outputs.initOutputs();
        end
        
        function obj = preprocessData(obj)
            % Preprocess video data.
            
            fprintf("\n----------------------------------\nVideo PreProcessing\n----------------------------------\n");
            
            obj.M0_data_video = obj.M0_raw_video;
            obj.M0_ff_video = obj.M0_ff_raw_video;
            obj.M1_data_video = obj.M1_raw_video;
            obj.M2_data_video = obj.M2_raw_video;
            
            % Register video
            tic;
            fprintf("\n----------------------------------\nVideo Registering\n----------------------------------\n");
            obj = VideoRegistering(obj);
            fprintf("- Video Registering took: %ds\n", round(toc));
            
            % Crop video
            tic;
            fprintf("\n----------------------------------\nVideo Cropping\n----------------------------------\n");
            obj = VideoCropping(obj);
            fprintf("- Video Cropping took: %ds\n", round(toc));
            
            % Normalize moments
            tic;
            fprintf("\n----------------------------------\nMoment Normalizing\n----------------------------------\n");
            obj = VideoNormalizingLocally(obj);
            fprintf("- Moment Normalizing took: %ds\n", round(toc));
            
            % Resize video
            tic;
            fprintf("\n----------------------------------\nVideo Resizing\n----------------------------------\n");
            obj = VideoResizing(obj);
            fprintf("- Video Resizing took: %ds\n", round(toc));
            
            % Interpolate video
            tic;
            fprintf("\n----------------------------------\nVideo Interpolation\n----------------------------------\n");
            obj = VideoInterpolating(obj);
            fprintf("- Video Interpolation took: %ds\n", round(toc));
            
            % Remove outliers
            tic;
            fprintf("\n----------------------------------\nVideo Outlier Cleaning\n----------------------------------\n");
            obj = VideoRemoveOutliers(obj);
            fprintf("- Video Outlier Cleaning took: %ds\n", round(toc));
            
            obj.is_preprocessed = true;

            obj.Outputs.initOutputs();

        end
        
        function obj = analyzeData(obj, app)
            % Main routine for EyeFlow analysis.
            
            % Initialize ToolBox and parameters
            ToolBox = obj.ToolBoxMaster;
            params = ToolBox.getParams;
            veins_analysis = params.veins_analysis;
            totalTime = tic;
            saveGit;
            ToolBox.Outputs = obj.Outputs;
            ToolBox.Outputs.add('NumFrames', size(obj.M0_data_video, 3), '', 0);
            ToolBox.Outputs.add('FrameRate', ToolBox.fs * 1000 / ToolBox.stride , 'Hz', 0);
            ToolBox.Outputs.add('InterFramePeriod', ToolBox.stride / ToolBox.fs / 1000 , 's', 0);
            
            % Mask Creation
            if obj.flag_segmentation
                fprintf("\n----------------------------------\nMask Creation\n----------------------------------\n");
                createMasksTimer = tic;
                
                f_AVG_mean = squeeze(mean(obj.f_AVG_video, 3));
                [obj.maskArtery, obj.maskVein, obj.maskNeighbors, obj.xy_barycenter] = ...
                    createMasks(obj.M0_ff_video, f_AVG_mean);
                
                M0_ff_img = rescale(mean(obj.M0_ff_video, 3));
                cmapArtery = ToolBox.cmapArtery;
                cmapVein = ToolBox.cmapVein;
                cmapAV = ToolBox.cmapAV;
                
                M0_Artery = setcmap(M0_ff_img, obj.maskArtery, cmapArtery);
                M0_Vein = setcmap(M0_ff_img, obj.maskVein, cmapVein);
                M0_AV = setcmap(M0_ff_img, obj.maskArtery & obj.maskVein, cmapAV);
                
                M0_RGB = (M0_Artery + M0_Vein) .* ~(obj.maskArtery & obj.maskVein) + M0_AV + rescale(M0_ff_img) .* ~(obj.maskArtery | obj.maskVein);
                app.ImageDisplay.ImageSource = mat2gray(M0_RGB); % Rescale the image for display
                ax = ancestor(app.ImageDisplay, 'axes');
                axis(ax, 'equal');
                
                obj.is_segmented = true;
                
                fprintf("- Mask Creation took: %ds\n", round(toc(createMasksTimer)));
            end
            
            % Pulse Analysis
            if obj.flag_bloodFlowVelocity_analysis
                fprintf("\n----------------------------------\nFind Systole\n----------------------------------\n");
                
                if ~isfolder(fullfile(ToolBox.path_png, 'bloodFlowVelocity'))
                    mkdir(ToolBox.path_png, 'bloodFlowVelocity')
                    mkdir(ToolBox.path_eps, 'bloodFlowVelocity')
                end
                
                fprintf("\n----------------------------------\nBlood Flow Velocity Analysis\n----------------------------------\n");
                pulseAnalysisTimer = tic;
                
                f_AVG_mean = squeeze(mean(obj.f_AVG_video, 3));
                [obj.vRMS, obj.sysIdxList, obj.sysIdx, obj.diasIdx] = pulseAnalysis(obj.f_RMS_video, obj.maskArtery, obj.maskVein, obj.maskNeighbors, obj.xy_barycenter);
                
                if params.json.PulseAnalysis.ExtendedFlag
                    extendedPulseAnalysis(obj.M0_ff_video, obj.f_RMS_video, f_AVG_mean, obj.vRMS, obj.maskArtery, obj.maskVein, obj.xy_barycenter, obj.sysIdxList);
                end
                
                obj.is_pulseAnalyzed = true;
                
                fprintf("- Blood Flow Velocity Analysis took: %ds\n", round(toc(pulseAnalysisTimer)));
            end
            
            % Pulse Velocity Analysis
            %  if obj.flag_pulseVelocity_analysis
            %     fprintf("\n----------------------------------\nPulse Velocity Calculation\n----------------------------------\n");
            %     pulseVelocityTimer = tic;
            
            %     pulseVelocity(obj.M0_data_video, maskArtery)
            
            %     time_pulsevelocity = toc(pulseVelocityTimer);
            %     fprintf("- Pulse Velocity Calculations took : %ds\n", round(time_pulsevelocity))
            % end
            
            % Blood Flow Velocity Analysis
            if obj.flag_bloodFlowVelocity_figures
                fprintf("\n----------------------------------\nBlood Flow Velocity Figures\n----------------------------------\n");
                bloodFlowVelocityTimer = tic;
                
                bloodFlowVelocity(obj.vRMS, obj.maskArtery, obj.maskVein, obj.M0_ff_video, obj.xy_barycenter);
                
                fprintf("- Blood Flow Velocity Figures calculation took: %ds\n", round(toc(bloodFlowVelocityTimer)));
            end
            
            % Cross-Section Analysis
            if obj.flag_crossSection_analysis
                fprintf("\n----------------------------------\nCross-Section Analysis\n----------------------------------\n");
                crossSectionAnalysisTimer = tic;
                
                [obj.Q_results_A] = crossSectionsAnalysis(obj.maskArtery, 'Artery', obj.vRMS, obj.M0_ff_video, obj.xy_barycenter);
                
                if veins_analysis
                    [obj.Q_results_V] = crossSectionsAnalysis(obj.maskVein, 'Vein', obj.vRMS, obj.M0_ff_video, obj.xy_barycenter);
                end
                
                obj.is_crossSectionAnalyzed = true;
                
                fprintf("- Cross-Section Analysis took: %ds\n", round(toc(crossSectionAnalysisTimer)));
            end
            
            % Cross-Section Figures
            if obj.flag_crossSection_figures
                fprintf("\n----------------------------------\nCross-Section Figures\n----------------------------------\n");
                crossSectionFiguresTimer = tic;
                
                crossSectionsFigures(obj.Q_results_A, obj.maskArtery, 'Artery', obj.M0_ff_video, obj.xy_barycenter, obj.sysIdxList, obj.sysIdx, obj.diasIdx);
                
                if veins_analysis
                    crossSectionsFigures(obj.Q_results_V, obj.maskVein, 'Vein', obj.M0_ff_video, obj.xy_barycenter, obj.sysIdxList, obj.sysIdx, obj.diasIdx);
                end
                
                try
                    generateHealthReport()
                catch ME
                    fprintf("Error generating health report: %s\n", ME.message);
                    for i = 1:length(ME.stack)
                        fprintf("Error in %s at line %d: %s\n", ME.stack(i).name, ME.stack(i).line, ME.message);
                    end
                end
                
                fprintf("- Cross-Section Figures took: %ds\n", round(toc(crossSectionFiguresTimer)));
            end
            
            % Spectral Analysis
            if obj.flag_spectral_analysis && isfile(fullfile(ToolBox.EF_path, 'raw', [strcat(ToolBox.main_foldername, '_SH'), '.raw']))
                fprintf("\n----------------------------------\nSpectral Analysis\n----------------------------------\n");
                timeSpectralAnalysis = tic;
                
                % Import SH data
                tmpname = strcat(ToolBox.main_foldername, '_SH');
                fileID = fopen(fullfile(ToolBox.EF_path, 'raw', [tmpname, '.raw']));
                videoSH = fread(fileID, 'float32');
                fclose(fileID);
                
                [numX, numY, numFrames] = size(obj.M0_data_video);
                bin_x = 4; bin_y = 4; bin_t = 1;
                SH_cube = reshape(videoSH, ceil(numX / (bin_x)), ceil(numY / (bin_y)), [], ceil(numFrames / bin_t));
                
                % Spectrum Analysis
                fprintf("\n----------------------------------\nSpectrum Analysis\n----------------------------------\n");
                spectrumAnalysisTimer = tic;
                
                spectrum_analysis(SH_cube, obj.M0_data_video);
                
                time_spectrumAnalysis = toc(spectrumAnalysisTimer);
                fprintf("- Spectrum Analysis took : %ds\n", round(time_spectrumAnalysis))
                
                % Spectrogram
                fprintf("\n----------------------------------\nSpectrogram\n----------------------------------\n");
                spectrogramTimer = tic;
                
                spectrum_video(SH_cube, obj.maskArtery, obj.maskNeighbors);
                
                fprintf("- Spectrogram took: %ds\n", round(toc(spectrogramTimer)));
                
                fprintf("\n----------------------------------\nSpectral Analysis timing: %ds\n", round(toc(timeSpectralAnalysis)));
            end
            
            % Main Outputs Saving
            
            fid = fopen(fullfile(ToolBox.path_json, strcat(ToolBox.main_foldername, '_EF_main_outputs.json')), 'w');
            fwrite(fid, jsonencode(ToolBox.outputs, "PrettyPrint", true), 'char');
            fclose(fid);
            
            ToolBox.Outputs.writeJson(fullfile(ToolBox.path_json, strcat(ToolBox.folder_name, '_main_outputs.json')));
            
            % Final Output
            tTotal = toc(totalTime);
            fprintf("\n----------------------------------\nTotal EyeFlow timing: %ds\n", round(tTotal));
            fprintf("End Computer Time: %s\n", datetime('now', 'Format', 'yyyy/MM/dd HH:mm:ss'));
            
            clear ToolBox;
            diary off;
            displaySuccessMsg(1);
        end
        
    end
    
end
