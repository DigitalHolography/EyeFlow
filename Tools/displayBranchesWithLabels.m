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

    for i = 1:N
        mask = (vesselBranchesMap==i);
        s = regionprops(mask,"centroid");
        locsLabels{i} = cat(1,s.Centroid);

        M0_RGB = M0_RGB .* mask .* cmap(i) + M0_RGB .* ~mask;
    end


    displayLabels(M0_RGB,locsLabels,1:N,'save_path',opt.save_path);

end