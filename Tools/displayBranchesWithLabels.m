function displayBranchesWithLabels(vesselBranchesMap,opt)

    arguments
        vesselBranchesMap
        opt.save_path = [];
        opt.bkgimg = [];
        opt.label = [] % 
        opt.unit = []
    end

    ToolBox = getGlobalToolBox();
    M0_RGB = ToolBox.Cache.M0_RGB;
    N = max(vesselBranchesMap(:)); % num of branches

    % cmap = colormap("hsv");
    

    figure("Visible","off")

    for i = 1:N
        mask = (vesselBranchesMap==i);
        s = regionprops(mask,"centroid");
        locsLabels{i} = cat(1,s.Centroid);
    end

    if ~isempty(opt.bkgimg)
        bkgimg = opt.bkgimg;
    else
        bkgimg = M0_RGB;
    end

    if isempty(opt.label)
        opt.label = 1:N;
    end


    displayLabels(bkgimg,locsLabels,opt.label,'save_path',opt.save_path,'unit',opt.unit);

end