function nonRigidMaskWithFrangiVesselness(refImg, targImg, maskIn, maskOut)
[~, frangi1] = frangiVesselness(refImg);
betterfrangi1 = (mat2gray(im2double(refImg) + 5000 * frangi1));

[~, frangi2] = frangiVesselness(targImg);
betterfrangi2 = (mat2gray(im2double(targImg) + 5000 * frangi2));

imshowpair(betterfrangi1, betterfrangi2);
nonRigidMask(betterfrangi1, betterfrangi2, maskIn, maskOut);
