function exportSegmentsValueToH5(name, maskLabel, values, baseName, units)
    arguments
        name 
        maskLabel 
        values 
        baseName string = "Segments"
        units string = ""
    end
    
    if ~endsWith(baseName, "/") && ~endsWith(baseName, "_")
        baseName = baseName + "/";
    end
    
    ToolBox = getGlobalToolBox;
    
    [numCircles, numBranches] = size(maskLabel);
    typeName = class(values(1));
    
    L = zeros(size(maskLabel{1,1}), typeName);
    % Process each circle and branch
    for cIdx = 1:numCircles
        for bIdx = 1:numBranches
            L(maskLabel{cIdx, bIdx}) = values(cIdx, bIdx);
        end
    end
    
    ToolBox.Output.add(sprintf("%s_Segments_Labels", name), L, h5path = baseName + sprintf("%s_Segments_Labels", name), unit = units);

end