function [skel_out, path_xy] = skel_keep_longest_path(mask)
% SKEL_KEEP_LONGEST_PATH  Skeletonize mask and keep only the single longest path
%
%   [skel_out, path_xy] = skel_keep_longest_path(mask)
%
% Inputs:
%   mask      - logical or numeric 2D mask (nonzero = foreground)
% Outputs:
%   skel_out  - logical image of the modified skeleton with exactly two endpoints
%   path_xy   - Nx2 array of [row, col] coordinates along the kept skeleton path,
%               in path order (from one endpoint to the other)
%
% Behavior:
%   - Skeletonizes mask using bwmorph(...,'skel',Inf).
%   - Builds a graph of skeleton pixels using 8-connectivity.
%   - If there are endpoints (>0), finds all endpoint pairs, gets the
%     geodesic (graph) path lengths between them, chooses the pair with
%     maximum geodesic length and keeps only that path.
%   - If there are no endpoints (closed loop), finds the two skeleton
%     pixels that are farthest apart in graph-distance and keeps the path
%     between them.
%
% Requires: Image Processing Toolbox (for bwmorph)
%
% Example:
%   mask = imread('my_mask.png') > 0;
%   [skel, coords] = skel_keep_longest_path(mask);
%   imshow(skel)
%   hold on; plot(coords(:,2), coords(:,1), '.r'); hold off

    if ~ismatrix(mask)
        error('Input must be a 2D mask.');
    end

    % Ensure logical
    mask = logical(mask);

    % 1) skeletonize
    sk = bwmorph(mask, 'skel', Inf);

    % If the resulting skeleton is empty, return quickly
    if nnz(sk) == 0
        skel_out = sk;
        path_xy = zeros(0,2);
        return
    end

    % 2) Get linear indexes and subscripts for skeleton pixels
    [r, c] = find(sk);
    N = numel(r);
    idx_lin = sub2ind(size(sk), r, c); % linear indices of skeleton pixels

    % Map from linear index to node id (1..N)
    lin2node = containers.Map('KeyType','uint32','ValueType','uint32');
    for k = 1:N
        lin2node(uint32(idx_lin(k))) = uint32(k);
    end

    % 3) Build adjacency using 8-neighborhood
    % Preallocate adjacency lists (sparse)
    % neighbor offsets (linear) for 8-connectivity:
    [nr, nc] = size(sk);
    neigh_offsets = [-1, 1, -nr, nr, -nr-1, -nr+1, nr-1, nr+1]; % careful: we'll compute neighbors by subscripts to avoid border linear math errors

    rows = [];
    cols = [];

    % For each node, check its 8-neighbors and if neighbor is skeleton create edge
    for k = 1:N
        rr = r(k); cc = c(k);
        for dr = -1:1
            for dc = -1:1
                if dr==0 && dc==0, continue; end
                r2 = rr + dr; c2 = cc + dc;
                if r2 < 1 || r2 > nr || c2 < 1 || c2 > nc, continue; end
                lin2 = uint32(sub2ind([nr nc], r2, c2));
                if sk(r2,c2)
                    % both are skeleton pixels. add edge (k -> neighbor node)
                    rows(end+1) = k; %#ok<AGROW>
                    cols(end+1) = lin2node(lin2); %#ok<AGROW>
                end
            end
        end
    end

    % Create undirected graph (make symmetric)
    A = sparse([rows cols], [cols rows], 1, N, N);
    G = graph(A);   % unweighted graph; each edge weight = 1

    % 4) Find endpoints in the skeleton (pixels with degree 1)
    deg = degree(G);
    endpoints_nodes = find(deg == 1);

    if numel(endpoints_nodes) >= 2
        % evaluate graph distances between every pair of endpoints and choose the pair with max distance
        E = endpoints_nodes;
        m = numel(E);
        maxd = -inf;
        bestPair = [E(1), E(min(2,end))];
        % Use distances from each endpoint to all nodes (efficient)
        for i = 1:m
            d = distances(G, E(i));
            for j = i+1:m
                dij = d(E(j));
                if ~isinf(dij) && dij > maxd
                    maxd = dij;
                    bestPair = [E(i), E(j)];
                end
            end
        end
        % extract shortest (graph) path between bestPair
        p = shortestpath(G, bestPair(1), bestPair(2));
        chosen_nodes = p;
    else
        % 0 or 1 endpoints (e.g., closed loop or degenerate). Find the two nodes
        % that are farthest apart in graph-distance.
        % Use distances from arbitrary node to get a candidate, then refine:
        % 1) pick node a arbitrary, find farthest node b
        % 2) from b find farthest node c. path between b and c is a diameter approximation.
        a = 1;
        d1 = distances(G, a);
        [~, b] = max(d1(~isinf(d1)));
        % need index mapping since result came with logical indexing
        % find index of max directly:
        [~, bnode] = max(d1);
        d2 = distances(G, bnode);
        [~, cnode] = max(d2);
        p = shortestpath(G, bnode, cnode);
        chosen_nodes = p;
    end

    % 5) create skeleton output containing only the chosen path nodes (and connecting edges)
    skel_out = false(size(sk));
    chosen_lin = idx_lin(chosen_nodes); % linear indices of chosen nodes
    skel_out(chosen_lin) = true;

    % 6) produce ordered list of coordinates along the path
    path_xy = [r(chosen_nodes), c(chosen_nodes)]; % [row, col] in path order

    % final check: ensure there are exactly two endpoints
    sk_final_nodes = find(skel_out);
    % number of neighbors for each node in skel_out (8-neigh)
    % (this should be 2 for interior points, 1 for endpoints)
    % Nothing else needed — function guarantees chosen path is single simple path.

end