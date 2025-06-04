function [maskLongArtery, label, adjMatrix] = getLongestArteryBranch(mask, xy_barycenter, name)
% Returns the longest artery branch mask
ToolBox = getGlobalToolBox;
params = ToolBox.getParams;
[numX, numY] = size(mask);

maskRadius = params.json.Mask.CropChoroidRadius;

x_barycenter = xy_barycenter(1);
y_barycenter = xy_barycenter(2);

%% Skeletonize and label the individual branches
skel = bwskel(mask);
[x, y] = meshgrid(1:numY, 1:numX);
cercleMask = sqrt((x - x_barycenter) .^ 2 + (y - y_barycenter) .^ 2) <= maskRadius / 6 * (numY + numX) / 2; % making a mask to cut the center

skel = skel & ~cercleMask; % takes out arteries near the center
[label, n] = bwlabel(skel & ~imdilate(bwmorph(skel, 'branchpoints'), strel('disk', 2))); % labeling the skeleton with branch points off to get individual branches

% FIX ME we need a graph based approach to get full successive branches.

% simply getting the longuest individual branch for now
mx = 0;
index = 1;

for i = 1:n % for each individual branch
    s = sum(label == i, 'all'); % number of pixel of this branch
    
    if s > mx
        mx = s;
        index = i;
    end
    
end

%% Get the label mask back to initial size
for i = 1:n % for each individual branch
    sk_mask = label == i;
    label(bwareafilt(imdilate(sk_mask, strel('disk', floor(numX * maskRadius / 10))) & mask, 1)) = i;
end

figure(Visible="off");
imagesc(label);
foldername = ToolBox.folder_name;
cmap = jet(max(label, [], 'all') + 1);
imwrite(label + 1, cmap, fullfile(ToolBox.path_png, 'mask', sprintf("%s_maskSeparateBranchesLabeled_%s.png", foldername, name)), 'png');

%% compute the adjency matrix of the artery tree and make the graph
adjMatrix = false(n);
skel2 = bwskel(mask);
bp = find(bwmorph(skel2, 'branchpoints'));

all_circles = false(size(mask)); %display only
xy_barycenter_coord = xy_barycenter(1) * numX + xy_barycenter(2);
circle = false(size(mask));
circle(xy_barycenter_coord) = true;
circle = imdilate(circle, strel('diamond', 8));
all_circles = all_circles | circle;

for i = 1:length(bp) % iterates over all branch points
    circle = false(size(mask));
    circle(bp(i)) = true;
    circle = imdilate(circle, strel('disk', 10));
    all_circles = all_circles | circle;
    neighboors = find(label & circle);
    
    for j = 1:length(neighboors) % iterates over all the neigbors of that branch point
        
        for k = 1:length(neighboors)
            adjMatrix(label(neighboors(j)), label(neighboors(k))) = true;
        end
        
    end
    
end

fi=figure(Visible="off"); imshow(all_circles | skel2); saveas(fi, fullfile(ToolBox.path_png, 'mask', sprintf("%s_maskTree_%s.png", foldername, name)));
g = graph(adjMatrix); fi=figure(Visible="off"); plot(g); saveas(fi, fullfile(ToolBox.path_png, 'mask', sprintf("%s_graphTree_%s.png", foldername, name)));
%% get the graph's biggest connected component
conn = conncomp(g);
idx_most_frequent_conn = mode(conn);

%% get the closest branch to the CRA
closest_branch = 0;
dist_closest_branch = Inf;

for i = 1:length(bp) % iterates over all branch points
    [xx, yy] = ind2sub(size(mask), bp(i));
    
    if (xx - x_barycenter) ^ 2 + (yy - y_barycenter) ^ 2 < dist_closest_branch && conn(label(bp(i))) == idx_most_frequent_conn
        dist_closest_branch = (xx - x_barycenter) ^ 2 + (yy - y_barycenter) ^ 2;
        closest_branch = label(bp(i));
    end
    
end

%% get the shortest path between branches from the closest branch to the CRA gathering the most pixels

% d = distances(g); % get the distance between each branch
%
% [~,i] = max(d(closest_branch,d(closest_branch,:)<Inf));
% P = shortestpath(g,closest_branch,i);
%
% maskLongArtery = false(size(maskArtery));
% for idx= 1:n
%     if conn(idx) == idx_most_frequent_conn && ismember(idx,P)
%         maskLongArtery = maskLongArtery | (L==idx);
%     end
% end
max_pixels_intensity = 0;
end_path = [];

for idx = 1:n
    
    if conn(idx) == idx_most_frequent_conn
        P = shortestpath(g, closest_branch, idx);
        pixels_intensity = 0;
        
        for j = 1:length(P)
            pixels_intensity = pixels_intensity + sum((label == P(j)), [1, 2]);
        end
        
        if pixels_intensity > max_pixels_intensity
            max_pixels_intensity = pixels_intensity;
            end_path = P;
        end
        
    end
    
end

maskLongArtery = false(size(mask));

for idx = 1:n
    
    if conn(idx) == idx_most_frequent_conn && ismember(idx, end_path)
        maskLongArtery = maskLongArtery | (label == idx);
    end
    
end

%%FIX ME going back to simplest option for reliability
m = 0;

for idx = 1:n
    
    if sum((label == idx), [1, 2]) > m
        maskLongArtery = (label == idx);
        m = sum((label == idx), [1, 2]);
    end
    
end

fi=figure(Visible="off"); imshow(maskLongArtery);

end
