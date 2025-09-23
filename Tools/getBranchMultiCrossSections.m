function [U_x,L] = getBranchMultiCrossSections(U,mask,absx,absy,halfwidth,numinterp,win)
% getBranchMultiCrossSections
% Computes cross-sectional profiles along a curve with smoothed tangents.
%
% Inputs:
%   U          - image stack (HxWxT)
%   mask       - binary mask (HxW)
%   absx, absy - coordinates along curve
%   halfwidth  - half width of strip
%   numinterp  - number of interpolation samples along strip
%   win        - smoothing window (default: 10 points)
%
% Outputs:
%   U_x - sampled intensity profiles (points x frames x numinterp)
%   L   - labeled strip mask

if nargin < 7
    win = 10; % default moving average window
end

numpoints = length(absx);
numFrames = size(U,3);

L   = zeros(size(mask));
U_x = zeros(numpoints, size(U,3), numinterp);

prev_line = [];

for i = 2:numpoints-1
    % --- compute averaged tangent ---
    imin = max(1, i-win);
    imax = min(numpoints, i+win);
    
    tx = absx(imin+1:imax) - absx(imin:imax-1);
    ty = absy(imin+1:imax) - absy(imin:imax-1);
    
    % remove degenerate steps
    valid = (tx.^2 + ty.^2) > 0;
    if ~any(valid), continue; end
    
    tvec = [mean(tx(valid)), mean(ty(valid))];
    if norm(tvec)==0, continue; end
    tangent = tvec / norm(tvec);
    normal  = [-tangent(2), tangent(1)];

    % endpoints of current line
    P3 = [absx(i) - halfwidth*normal(1), absy(i) - halfwidth*normal(2)];
    P4 = [absx(i) + halfwidth*normal(1), absy(i) + halfwidth*normal(2)];

    % --- build quadrilateral strip between consecutive lines ---
    if ~isempty(prev_line)
        P1 = prev_line(1,:); % previous start
        P2 = prev_line(2,:); % previous end

        % parameter grid (resolution of strip filling)
        nu = 4*(2*halfwidth+1);   % 4x oversampling
        nv = 4*round(sqrt(sum((P3-P1).^2)));
        [u,v] = meshgrid(linspace(0,1,nu), linspace(0,1,nv));

        % bilinear interpolation of quadrilateral
        X = (1-v).*((1-u)*P1(1) + u*P2(1)) + v.*((1-u)*P3(1) + u*P4(1));
        Y = (1-v).*((1-u)*P1(2) + u*P2(2)) + v.*((1-u)*P3(2) + u*P4(2));

        % round to pixel indices
        X = round(X); 
        Y = round(Y);

        % keep inside image/mask
        inside = X>=1 & X<=size(mask,2) & Y>=1 & Y<=size(mask,1);
        X = X(inside); 
        Y = Y(inside);
        idx = sub2ind(size(mask), Y, X);

        if ~isempty(idx)
            L(idx) = i;
            for j=1:numFrames
                idx_t = sub2ind(size(U), Y(:), X(:), repelem(j, length(Y))');
                profile = U(idx_t);
                % normalize x-axis to [0,1]
                xi = (1:length(profile))/length(profile);
                U_x(i,j,:) = interp1(xi, profile, linspace(0,1,numinterp), 'linear', 0);
            end
        end
    end

    prev_line = [P3; P4]; % save endpoints for next step
end
end
