function pulseVelocity(M, maskArtery)
ToolBox = getGlobalToolBox;
[maskLongArtery, L, adjMatrix] = getLongestArteryBranch(maskArtery, ToolBox.Cache.list.xy_barycenter,'artery');

%maskLongArtery = (L==10);
PWV = pulseWaveVelocity(M,maskLongArtery);
end
