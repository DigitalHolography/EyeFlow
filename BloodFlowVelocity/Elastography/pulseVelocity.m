function pulseVelocity(M, ~, maskVessel, name)
ToolBox = getGlobalToolBox;

[L, n] = labelVesselBranches(maskVessel, ones(size(maskVessel)), ToolBox.Cache.xy_barycenter);

PWV = NaN(1, n);

for i = 1:n
    % displacementAnalysis(D, maskLongArtery);
    PWV(i) = pulseWaveVelocity(M, L == i, i, name);
end

close all;
end
