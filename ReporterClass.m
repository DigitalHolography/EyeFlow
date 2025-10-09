classdef ReporterClass < handle
% Handles output generation and reporting

properties

end

methods

    function obj = ReporterClass(executionObj)
        toolbox = getGlobalToolBox;

        if isempty(toolbox)
            error("ToolBoxMaster is not initialized in ExecutionClass.");
        end

        output = executionObj.Output;

        if isempty(output)
            error("Output is not initialized in ExecutionClass.");
        end

        toolbox.Output.add('NumFrames', size(executionObj.M0, 3), '', 0);
        toolbox.Output.add('FrameRate', toolbox.fs * 1000 / toolbox.stride, 'Hz', 0);
        toolbox.Output.add('InterFramePeriod', toolbox.stride / toolbox.fs / 1000, 's', 0);

        if ~isempty(toolbox.record_time_stamps_us)
            tmp = toolbox.record_time_stamps_us;
            toolbox.Output.add('UnixTimestampFirst', tmp.first, 'µs');
            toolbox.Output.add('UnixTimestampLast', tmp.last, 'µs');
        end

        if ~isfile(fullfile(toolbox.path_gif, sprintf("%s_M0.gif", toolbox.folder_name)))
            writeGifOnDisc(imresize(rescale(executionObj.M0_ff), 0.5), "M0")
        end

    end

    function saveOutputs(~, ~)
        ToolBox = getGlobalToolBox;
        ToolBox.Output.writeJson(fullfile(ToolBox.path_json, strcat(ToolBox.folder_name, 'output.json')));
        ToolBox.Output.writeHdf5(fullfile(ToolBox.path_h5, strcat(ToolBox.folder_name, 'output.h5')));
    end

    function getA4Report(~, executionObj)

        ToolBox = getGlobalToolBox;

        if executionObj.is_pulseAnalyzed && executionObj.veins_analysis

            try
                generateA4Report();
            catch ME
                MEdisp(ME, ToolBox.EF_path);
            end

        end

    end

    function displayFinalSummary(~, totalTime)
        tTotal = toc(totalTime);
        fprintf("\n----------------------------------\nTotal EyeFlow timing: %ds\n", round(tTotal));
        fprintf("End Computer Time: %s\n", datetime('now', 'Format', 'yyyy/MM/dd HH:mm:ss'));
        displaySuccessMsg(1);
    end

end

end
