function perBeatAnalysis()
    ToolBox = getGlobalToolBox;
    params = ToolBox.params;
    dt = ToolBox.stride / ToolBox.fs / 1000;

    bandLimitedSignalHarmonicCount  = params.json.PulseAnalysis.BandLimitedSignalHarmonicCount;
    sys_idx_list                    = ToolBox.Cache.sysIdxList;
    v_vein_signal                   = ToolBox.Cache.VeinVelocity;
    v_artery_signal                 = ToolBox.Cache.ArterialVelocity;

    ToolBox.Output.add("beatPeriodIdx", diff(sys_idx_list), h5path = "/Artery/PerBeat/beatPeriodIdx");
    ToolBox.Output.add("beatPeriodSeconds", diff(sys_idx_list) * dt, h5path = "/Artery/PerBeat/beatPeriodSeconds");

    perBeatAnalysis_handle(v_vein_signal,   "vein",     sys_idx_list, bandLimitedSignalHarmonicCount);
    perBeatAnalysis_handle(v_artery_signal, "artery",   sys_idx_list, bandLimitedSignalHarmonicCount);
end

function perBeatAnalysis_handle(v_signal, name, sys_idx_list, bandLimitedSignalHarmonicCount)
    arguments
        v_signal
        name string
        sys_idx_list
        bandLimitedSignalHarmonicCount
    end

    ToolBox = getGlobalToolBox;
    dt = ToolBox.stride / ToolBox.fs / 1000;
        
    [VelocitySignalPerBeat, VelocitySignalPerBeatFFT, VelocitySignalPerBeatBandLimited] = perBeatSignalAnalysis(v_signal, sys_idx_list, bandLimitedSignalHarmonicCount);

    ToolBox.Output.add(name + "_VelocitySignalPerBeat",             VelocitySignalPerBeat,               h5path = capitalize(name) + "/PerBeat/VelocitySignalPerBeat",              keepSize=true);
    ToolBox.Output.add(name + "_VelocitySignalPerBeatFFT_abs",      abs(VelocitySignalPerBeatFFT),       h5path = capitalize(name) + "/PerBeat/VelocitySignalPerBeatFFT_abs",       keepSize=true);
    ToolBox.Output.add(name + "_VelocitySignalPerBeatFFT_arg",      angle(VelocitySignalPerBeatFFT),     h5path = capitalize(name) + "/PerBeat/VelocitySignalPerBeatFFT_arg",       keepSize=true);
    ToolBox.Output.add(name + "_VelocitySignalPerBeatBandLimited",  VelocitySignalPerBeatBandLimited,    h5path = capitalize(name) + "/PerBeat/VelocitySignalPerBeatBandLimited",   keepSize=true);

    ToolBox.Output.add(name + "_VmaxPerBeatBandLimited", max(VelocitySignalPerBeatBandLimited, [], 2), h5path = capitalize(name) + "/PerBeat/VmaxPerBeatBandLimited");
    ToolBox.Output.add(name + "_VminPerBeatBandLimited", min(VelocitySignalPerBeatBandLimited, [], 2), h5path = capitalize(name) + "/PerBeat/VminPerBeatBandLimited");

    ToolBox.Output.add(name + "_VTIPerBeat", sum(VelocitySignalPerBeat, 2) * dt, h5path = capitalize(name) + "/PerBeat/VTIPerBeat");
end