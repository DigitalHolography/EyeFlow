function pulseVelocity(M,D, maskVessel,name)
ToolBox = getGlobalToolBox;

[L, n] = labelVesselBranches(maskVessel, ones(size(maskVessel)), ToolBox.Cache.list.xy_barycenter);

for i=1:n
    % displacementAnalysis(D, maskLongArtery);
    PWV(i) = pulseWaveVelocity(M,L==i,i,name);
end

end
