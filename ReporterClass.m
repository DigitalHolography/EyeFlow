classdef ReporterClass < handle
% Handles output generation and reporting

properties (Access = private)
    ToolBoxMaster
    Output
end

methods

    function obj = ReporterClass(toolbox, output)
        obj.ToolBoxMaster = toolbox;
        obj.Output = output;

        obj.Output.add('NumFrames', size(executionObj.M0, 3), '', 0);
        obj.Output.add('FrameRate', toolbox.fs * 1000 / toolbox.stride, 'Hz', 0);
        obj.Output.add('InterFramePeriod', toolbox.stride / toolbox.fs / 1000, 's', 0);

        if ~isempty(toolbox.record_time_stamps_us)
            tmp = toolbox.record_time_stamps_us;
            obj.Output.add('UnixTimestampFirst', tmp.first, 'µs');
            obj.Output.add('UnixTimestampLast', tmp.last, 'µs');
        end

        if ~isfile(fullfile(ToolBox.path_gif, sprintf("%s_M0.gif", ToolBox.folder_name)))
            writeGifOnDisc(imresize(rescale(executionObj.M0_ff), 0.5), "M0")
        end

    end

    function saveOutputs(obj)
        ToolBox = obj.ToolBoxMaster;
        obj.Output.writeJson(fullfile(ToolBox.path_json, strcat(ToolBox.folder_name, 'output.json')));
        obj.Output.writeHdf5(fullfile(ToolBox.path_h5, strcat(ToolBox.folder_name, 'output.h5')));
    end

    function generateA4Report(obj, executionObj)

        if executionObj.is_pulseAnalyzed && executionObj.veins_analysis

            try
                generateA4Report()
            catch ME
                MEdisp(ME, obj.ToolBoxMaster.EF_path)
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
