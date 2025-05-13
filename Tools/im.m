function im(A)
ToolBox = getGlobalToolBox;

if size(A, 3)
    A = mean(A, 3);
end

figure, imshow(A, [], Colormap = ToolBox.cmapArtery)
end
