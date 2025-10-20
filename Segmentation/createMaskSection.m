function [maskSection, VesselImageRGB] = createMaskSection(ToolBox, img, r1, r2, xy_barycenter, figname, maskArtery, maskVein, opt)

arguments
    ToolBox
    img
    r1
    r2
    xy_barycenter
    figname
    maskArtery = []
    maskVein = []
    opt.thin = 0
end

params = ToolBox.getParams;

if isempty(maskArtery) && isempty(maskVein)
    error('Add at least one mask')
end

if isempty(maskVein)
    [numX, numY] = size(maskArtery);
else
    [numX, numY] = size(maskVein);
end

x_c = xy_barycenter(1) / numX;
y_c = xy_barycenter(2) / numY;

if opt.thin > 0
    maskSection1 = diskMask(numX, numY, r1, r1 + opt.thin, center = [x_c, y_c]);
    maskSection2 = diskMask(numX, numY, r2 - opt.thin, r2, center = [x_c, y_c]);
    maskSection = xor(maskSection1, maskSection2);
else
    maskSection = diskMask(numX, numY, r1, r2, center = [x_c, y_c]);
end

VesselImageRGB = zeros(numX, numY, 3);

if params.json.save_figures
    % Create Colormap Artery/Vein
    img = mat2gray(img);
    img_RGB = cat(3, img, img, img);
    maskSectionRGB = cat(3, maskSection, maskSection, maskSection);

    if isempty(maskVein)

        maskArterySection = maskArtery & maskSection;

        cmapArtery = ToolBox.Cache.cmapArtery;
        cmapArterySection = flipud(cmapArtery);

        M0_Artery = setcmap(img, maskArtery, cmapArtery);
        M0_ArterySection = setcmap(img, maskArterySection, cmapArterySection);

        VesselImageRGB = M0_Artery .* ~maskSection + M0_ArterySection + ...
            maskSectionRGB .* ~maskArtery + ...
            img_RGB .* ~maskArtery .* ~maskSection;

        if ~isfolder(fullfile(ToolBox.path_png, 'mask'))
            mkdir(fullfile(ToolBox.path_png, 'mask'));
        end

        imwrite(VesselImageRGB, fullfile(ToolBox.path_png, 'mask', sprintf("%s_%s.png", ToolBox.folder_name, figname)), 'png');

    elseif isempty(maskArtery)

        maskVeinSection = maskVein & maskSection;

        cmapVein = ToolBox.Cache.cmapVein;
        cmapVeinSection = flipud(cmapVein);

        M0_Vein = setcmap(img, maskVein, cmapVein);
        M0_VeinSection = setcmap(img, maskVeinSection, cmapVeinSection);

        VesselImageRGB = M0_Vein .* ~maskSection + M0_VeinSection + ...
            maskSectionRGB .* ~maskVein + ...
            img_RGB .* ~maskVein .* ~maskSection;

        if ~isfolder(fullfile(ToolBox.path_png, 'mask'))
            mkdir(fullfile(ToolBox.path_png, 'mask'));
        end

        imwrite(VesselImageRGB, fullfile(ToolBox.path_png, 'mask', sprintf("%s_%s.png", ToolBox.folder_name, figname)), 'png');

    else

        maskArterySection = maskArtery & maskSection;
        maskVeinSection = maskVein & maskSection;
        maskAV = maskArtery & maskVein;
        maskAVSection = maskAV & maskSection;

        cmapArtery = ToolBox.Cache.cmapArtery;
        cmapArterySection = flipud(cmapArtery);
        cmapVein = ToolBox.Cache.cmapVein;
        cmapVeinSection = flipud(cmapVein);
        cmapAV = ToolBox.Cache.cmapAV;
        cmapAVSection = flipud(cmapAV);

        M0_Artery = setcmap(img, maskArtery, cmapArtery);
        M0_ArterySection = setcmap(img, maskArterySection, cmapArterySection);
        M0_Vein = setcmap(img, maskVein, cmapVein);
        M0_VeinSection = setcmap(img, maskVeinSection, cmapVeinSection);
        M0_AV = setcmap(img, maskAV, cmapAV);
        M0_AVSection = setcmap(img, maskAVSection, cmapAVSection);

        maskVessel = maskArtery | maskVein;

        VesselImageRGB = ...
            M0_Artery .* ~maskAV .* ~maskSection + ...
            M0_ArterySection .* ~maskAVSection + ...
            M0_Vein .* ~maskAV .* ~maskSection + ...
            M0_VeinSection .* ~maskAVSection + ...
            M0_AV .* ~maskSection + ...
            M0_AVSection + ...
            maskSectionRGB .* ~maskVessel + ...
            img_RGB .* ~maskArtery .* ~maskVein .* ~maskSection;

        if ~isfolder(fullfile(ToolBox.path_png, 'mask'))
            mkdir(fullfile(ToolBox.path_png, 'mask'));
        end

        imwrite(VesselImageRGB, fullfile(ToolBox.path_png, 'mask', sprintf("%s_%s.png", ToolBox.folder_name, figname)), 'png');
    end

end

end
