function displayBranchesWithLabels(vesselBranchesMap,opt)

    arguments
        vesselBranchesMap
        opt.save_path
    end

    ToolBox = getGlobalToolBox();
    M0_RGB = ToolBox.Cache.M0_RGB;
    N = max(vesselBranchesMap(:)); % num of branches

    cmap = colormap("hsv");
    

    figure("Visible","off")

    for i = 1:n
        mask = (vesselBranchesMap==i);
        s = regionprops(mask,"centroid");
        locsLabel{i} = cat(1,s.Centroid);

        M0_RGB = M0_RGB .* mask .* cmap(i) + M0_RGB .* ~mask;
    end


    displayLabels(M0_RGB,locsLabels,linspace(1,N),save_png_file_name=sprintf("%s_%s_branches"),opt.save_path);

end