function capturedFrame = displayLabels(BkgImg,locsLabels,Labels,opt)

    arguments
        BkgImg
        locsLabels
        Labels
        opt.save_path = []
        opt.unit = []
    end

    



    assert(numel(locsLabels)==numel(Labels));

    if ~iscell(Labels)
        Labels = num2cell(Labels);
    end

    Nlab = numel(Labels);
    fig = figure("Visible","off");
    imshow(BkgImg);
    fontsize = round(size(BkgImg,1) / 40);
    
    textOptions = {"FontWeight", "bold", ...
        "Color", "white", ...
        "FontSize", fontsize, ...
        "BackgroundColor", "black"};

    for i = 1:Nlab
        x_l = locsLabels{i}(1);
        y_l = locsLabels{i}(2);
        label = Labels{i};
        
        if isnumeric(label) && ~isnan(label)
            if ~isempty(opt.unit)
                text(x_l, y_l, sprintf("%.4g %s", round(label, 2), opt.unit), textOptions{:});
            else
                text(x_l, y_l, sprintf("%d", round(label, 1), opt.unit), textOptions{:});
            end
        elseif isstring(label)
            text(x_l, y_l, sprintf("%s", label), textOptions{:});
        end
    end

    if ~isempty(opt.save_path)
        ToolBox = getGlobalToolBox();
        saveas(gcf,opt.save_path);
    end


    capturedFrame = frame2im(getframe(gca));
end