function obj = readHDF5(obj)

dir_path_raw = fullfile(obj.directory, 'raw');
NameRefRawFile = strcat(obj.filenames, '_raw.h5');
RefRawFilePath = fullfile(dir_path_raw, NameRefRawFile);


try
    %h5disp(RefRawFilePath)

    

    obj.M0_raw_video = squeeze(h5read(RefRawFilePath,'/moment0'));
    obj.M1_raw_video = squeeze(h5read(RefRawFilePath,'/moment1'));
    obj.M2_raw_video = squeeze(h5read(RefRawFilePath,'/moment2'));
    %obj.M0_ff_raw_video = obj.M0_raw_video;


catch ME

    disp(['ID: ' ME.identifier])
    rethrow(ME)

end

try
    % Import SH data
    obj.SH_data_hypervideo = h5read(RefRawFilePath,'/SH');

catch
    warning('no SH was found')
end


end
