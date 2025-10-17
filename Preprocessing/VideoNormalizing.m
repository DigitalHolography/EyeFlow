function VideoNormalizing(obj)

params = Parameters_json(obj.directory, obj.param_name);
gwRatio = params.json.FlatFieldCorrection.GWRatio;
border = params.json.FlatFieldCorrection.Border;

M0_data_mean = mean(double(obj.M0), [1 2]);

obj.f_RMS = sqrt(double(obj.M2) ./ M0_data_mean);
obj.f_AVG = double(obj.M1) ./ M0_data_mean;
obj.M0_ff = flat_field_correction(obj.M0_ff, ceil(gwRatio * size(obj.M0_ff, 1)), border);

end
