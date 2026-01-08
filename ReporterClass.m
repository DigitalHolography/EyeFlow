classdef ReporterClass < handle
% Handles output generation and reporting

properties

end

methods

    function obj = ReporterClass(executionObj)
        ToolBox = getGlobalToolBox;
        params = ToolBox.getParams;

        if isempty(ToolBox)
            error("ToolBoxMaster is not initialized in ExecutionClass.");
        end

        output = executionObj.Output;

        if isempty(output)
            error("Output is not initialized in ExecutionClass.");
        end

        if isempty(executionObj.M0)
            ToolBox.Output.add('NumFrames', size(executionObj.M0_ff, 3), h5path = '/Meta/NumFrames');
        else
            ToolBox.Output.add('NumFrames', size(executionObj.M0, 3), h5path = '/Meta/NumFrames');
        end

        ToolBox.Output.add('FrameRate', ToolBox.fs * 1000 / ToolBox.stride, 'Hz', h5path = '/Meta/FrameRate');
        ToolBox.Output.add('InterFramePeriod', ToolBox.stride / ToolBox.fs / 1000, 's', h5path = '/Meta/InterFramePeriod');

        if ~isempty(ToolBox.record_time_stamps_us)
            tmp = ToolBox.record_time_stamps_us;
            ToolBox.Output.add('UnixTimestampFirst', tmp.first, 'µs', h5path = '/Meta/UnixTimestampFirst');
            ToolBox.Output.add('UnixTimestampLast', tmp.last, 'µs', h5path = '/Meta/UnixTimestampLast');
        end

        if ~isfile(fullfile(ToolBox.path_gif, sprintf("%s_M0.gif", ToolBox.folder_name))) && params.saveFigures
            writeGifOnDisc(imresize(rescale(executionObj.M0_ff), 0.5), "M0")
        end

    end

    function saveOutputs(~, ~)
        tic
        fprintf("Saving Outputs...\n");

        ToolBox = getGlobalToolBox;

        if ~isempty(ToolBox.Output)
            ToolBox.Output.writeJson(fullfile(ToolBox.path_json, sprintf("%s_output.json", ToolBox.folder_name)));
            disp("JSON output is done");
            ToolBox.Output.writeHdf5(fullfile(ToolBox.path_h5, sprintf("%s_output.h5", ToolBox.folder_name)));
            disp("H5 output is done");
        else
            DummyOutput = OutputClass();
            DummyOutput.writeJson(fullfile(ToolBox.path_json, sprintf("%s_output.json", ToolBox.folder_name)));
            DummyOutput.writeHdf5(fullfile(ToolBox.path_h5, sprintf("%s_output.h5", ToolBox.folder_name)));
        end

        fprintf("Saving Outputs took %ds\n", round(toc));
    end

    function getA4Report(~, ME)

        if nargin < 2
            ME = [];
        end

        if ~isempty(ME)
            ToolBox = getGlobalToolBox;
            str = MEdisp(ME, ToolBox.EF_path);

            fid = fopen(fullfile(ToolBox.path_log, sprintf("%s_error_log.txt", ToolBox.folder_name)), 'w');
            fprintf(fid, "%s", str);
            fclose(fid);
        end

        tic
        fprintf("Generating A4 Report...\n");

        ToolBox = getGlobalToolBox;

        if ~isempty(ToolBox.Output)
            ToolBox.Output.writeJson(fullfile(ToolBox.path_json, sprintf("%s_output.json", ToolBox.folder_name)));
            % ToolBox.Output.writeHdf5(fullfile(ToolBox.path_h5, sprintf("%s_output.h5", ToolBox.folder_name)));
        else
            ToolBox.Output = OutputClass();
        end

        try
            generateA4Report(ME);
        catch ME
            MEdisp(ME, ToolBox.EF_path);
        end

        fprintf("Generating A4 Report took %ds\n", round(toc));
    end

    function displayFinalSummary(~, totalTime)
        tTotal = toc(totalTime);
        fprintf("\n----------------------------------\nTotal EyeFlow timing: %ds\n", round(tTotal));
        fprintf("End Computer Time: %s\n", datetime('now', 'Format', 'yyyy/MM/dd HH:mm:ss'));
        displaySuccessMsg(1);
    end

    function generateReport(obj, executionObj)
        % Generate reports and outputs
        totalTime = tic;
        obj.getA4Report(executionObj);
        obj.saveOutputs();
        obj.displayFinalSummary(totalTime);

    end

end

end
