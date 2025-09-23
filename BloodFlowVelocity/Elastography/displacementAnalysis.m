function displacementAnalysis(displacementField, mask)

ToolBox = getGlobalToolBox;
params = ToolBox.getParams;

[numX, numY] = size(mask);

D = squeeze(complex(displacementField.field(:,:,1,:),displacementField.field(:,:,2,:)));
D = imresize(D,[numX,numY]);
numFrames = size(D, 3);

%%
[skel_out, path_xy] = skel_keep_longest_path(mask);

[U_x,L] = getBranchMultiCrossSections(D,skel_out,path_xy(:,1),path_xy(:,2),10,30,1);

figure(74);
imagesc(L)
title('selected sections along the artery')

numpoints = size(U_x,1);
numFrames = size(U_x,2);
Ux = zeros(numpoints, numFrames);

for k=1:numFrames
    for j=2:numpoints

        [r,lags]=xcorr(squeeze(U_x(j-1,k,:)),squeeze(U_x(j,k,:)),'unbiased');
        [~,ind] = max(r);
        % disp(lags(ind))
        Ux(j,k) = lags(ind);
    end
end

end