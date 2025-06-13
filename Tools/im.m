function im(A)
ToolBox = getGlobalToolBox;

if size(A, 3)
    A = mean(A, 3);
end

figure, imagesc(A),
axis image
colormap(ToolBox.cmapArtery)
colorbar
end
