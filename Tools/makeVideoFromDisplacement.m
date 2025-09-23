function mask_vid = makeVideoFromDisplacement(mask,D)

numFrames=size(D,4);

mask = imresize(mask,[size(D,1),size(D,2)]);

mask_vid = zeros([size(mask) numFrames]);


for i=1:numFrames
    mask_vid(:,:,i) = imwarp(mask, D(:,:,:,i), "nearest");
end

end
