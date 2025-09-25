function PWV = pulseWaveVelocity(U, mask)
% Computes the pulse wave velocity based on a cross correlation computation
% U is the field over which we compute the velocity and mask is the mask of
% the selected retinal artery

% U(x,y,t) usually M0
% center the [x,y] barycenter (the center of the CRA)
ToolBox = getGlobalToolBox;
params = ToolBox.getParams;

% mask = imerode(mask,strel('disk', 4));

% implay(rescale(abs(U)) .* mask);

[numX, numY] = size(mask);
N_frame = size(U, 3);

%% create a grid of points to select points along the skeletton of the artery mask
%
dxx = 5;
skel = bwskel(mask);
grid = ones([numY, numX]) < 0;
grid(1:dxx:end, :) = true;
grid(:, 1:dxx:end) = true;
interpoints = grid & skel; % get points interpolating with grid

numpoints = sum(interpoints, 'all'); % get the number of points

%% register the positions of points stating by the one closest to the CRA then going from closest to closest

[interpoints_y, interpoints_x] = ind2sub(size(interpoints), find(interpoints)); % y first
k = dsearchn([interpoints_x, interpoints_y], [x_bary, y_bary]); % get the nearest point to the center

absx = zeros([1, numpoints]); % x and y position register
absy = zeros([1, numpoints]);
abs_dist = zeros([1, numpoints]); % vessel curvilign absis

absx(1) = interpoints_x(k); % nearest point to the center
absy(1) = interpoints_y(k);
interpoints_x(k) = [];
interpoints_y(k) = [];
abs_dist(1) = 0;

for kb = 2:numpoints
    k = dsearchn([interpoints_x, interpoints_y], [absx(kb - 1), absy(kb - 1)]);
    absx(kb) = interpoints_x(k);
    absy(kb) = interpoints_y(k);
    interpoints_x(k) = []; % deleting point from list
    interpoints_y(k) = [];
end

for kb = 2:numpoints
    abs_dist(kb) = abs_dist(kb - 1) + sqrt((absx(kb) - absx(kb - 1)) ^ 2 + (absy(kb) - absy(kb - 1)) ^ 2) * params.px_size;
end


PWV = NaN;
%% create a grid of points to select points along the skeletton of the artery mask
%
dxx = 5;
skel = bwskel(mask);
grid = ones([numY, numX]) < 0;
grid(1:dxx:end, :) = true;
grid(:, 1:dxx:end) = true;
interpoints = grid & skel; % get points interpolating with grid

numpoints = sum(interpoints, 'all'); % get the number of points

%% register the positions of points stating by the one closest to the CRA then going from closest to closest

[interpoints_y, interpoints_x] = ind2sub(size(interpoints), find(interpoints)); % y first
k = dsearchn([interpoints_x, interpoints_y], [x_bary, y_bary]); % get the nearest point to the center

absx = zeros([1, numpoints]); % x and y position register
absy = zeros([1, numpoints]);
abs_dist = zeros([1, numpoints]); % vessel curvilign absis

absx(1) = interpoints_x(k); % nearest point to the center
absy(1) = interpoints_y(k);
interpoints_x(k) = [];
interpoints_y(k) = [];
abs_dist(1) = 0;

for kb = 2:numpoints
    k = dsearchn([interpoints_x, interpoints_y], [absx(kb - 1), absy(kb - 1)]);
    absx(kb) = interpoints_x(k);
    absy(kb) = interpoints_y(k);
    interpoints_x(k) = []; % deleting point from list
    interpoints_y(k) = [];
end

for kb = 2:numpoints
    abs_dist(kb) = abs_dist(kb - 1) + sqrt((absx(kb) - absx(kb - 1)) ^ 2 + (absy(kb) - absy(kb - 1)) ^ 2) * params.px_size;
end



%% also possible to do like this (chatgpt)
% Skeletonize
% sk = bwmorph(mask, 'skel', Inf);
% [y, x] = find(sk);
% 
% % Order skeleton points (approximate ordering using bwtraceboundary)
% B = bwtraceboundary(sk, [y(1) x(1)], 'N');
% absx = B(:,2); absy = B(:,1);
% numpoints = numel(absx);
% 
% % Compute arc length along skeleton
% dx = diff(absx); dy = diff(absy);
% abs_dist = [0; cumsum(sqrt(dx.^2 + dy.^2))]' * params.px_size;
% 
% figure(73)
% plot(abs_dist);
% Ltot = abs_dist(end);



%% for the positions extract a signal

% First idea with strel and imdilate but not really orthogonal
%
% st_el = strel('disk', floor(numX * params.json.Mask.DiaphragmRadius / 5));
% for i = 1:numpoints
%     sk_mask = false(size(mask));
%     sk_mask(absy(i), absx(i)) = true;
%     sectio = imdilate(sk_mask, st_el) & mask;
%     L(sectio) = i;
%     U_x(i, :) = squeeze(mean(U .* sectio, [1, 2]));
% end


% Second idea with orhtogonal sections
halfwidth = 10;

L   = zeros(size(mask));
numinterp=30;
U_x = zeros(numpoints, size(U,3),numinterp);

prev_line = [];

for i = 2:numpoints-1
    % tangent/normal
    tx = absx(i+1) - absx(i-1);
    ty = absy(i+1) - absy(i-1);
    if tx==0 && ty==0, continue; end
    tangent = [tx, ty] / norm([tx, ty]);
    normal  = [-tangent(2), tangent(1)];

    % endpoints of current line
    P3 = [absx(i) - halfwidth*normal(1), absy(i) - halfwidth*normal(2)];
    P4 = [absx(i) + halfwidth*normal(1), absy(i) + halfwidth*normal(2)];

    if ~isempty(prev_line)
        P1 = prev_line(1,:); % previous start
        P2 = prev_line(2,:); % previous end

        % parameter grid (controls resolution of strip filling)
        nu = 2*halfwidth+1;
        nv = round(sqrt(sum((P3-P1).^2))); % distance between strips
        [u,v] = meshgrid(linspace(0,1,nu), linspace(0,1,nv));

        % bilinear interpolation of quadrilateral
        X = (1-v).*((1-u)*P1(1) + u*P2(1)) + v.*((1-u)*P3(1) + u*P4(1));
        Y = (1-v).*((1-u)*P1(2) + u*P2(2)) + v.*((1-u)*P3(2) + u*P4(2));

        % round to pixel indices
        X = round(X); Y = round(Y);

        % keep inside image/mask
        inside = X>=1 & X<=size(mask,2) & Y>=1 & Y<=size(mask,1);
        X = X(inside); Y = Y(inside);
        idx = sub2ind(size(mask), Y, X);
        idx = idx(mask(idx));
        

        if ~isempty(idx)
            L(idx) = i;
            for j=1:N_frame
                idx_t = sub2ind(size(U), Y(:), X(:), repelem(j, length(Y))');
                U_x(i,j,:) = interp1((1:length(U(idx_t)))/length(U(idx_t)),U(idx_t),(1:numinterp)/numinterp);
            end
        end
    end

    prev_line = [P3; P4]; % save endpoints for next step
end

figure(74);
imagesc(L)
title('selected sections along the artery')

Ux = zeros(numpoints, size(U,3));

for k=1:N_frame
    for j=2:numpoints

        [r,lags]=xcorr(squeeze(U_x(j-1,k,:)),squeeze(U_x(j,k,:)),'unbiased');
        [~,ind] = max(r);
        % disp(lags(ind))
        Ux(j,k) = lags(ind);
    end
end

figure(75);
imagesc((Ux))

% ft_Ux = fft(Ux, [], 2);
% ph = angle(ft_Ux);

% figure(78);
% imagesc(log10(abs(ft_Ux)));
% figure(79);
% imagesc(ph);

% Ux = Ux - mean(Ux, 2);
% hUx = hilbert(Ux')';
Uy=mean(U_x,3);% exp(1j * angle(hUx'));
hUy = hilbert(Uy);
C = exp(1j * angle(hUy'));
R=xcorr(C,'unbiased');
figure(544), imagesc(real(C));
[Nlags, cols] = size(R);
M = size(C,2);
Ravg = zeros(Nlags, 2*M+1);
for i = -M:M
    idx = find_indices_compact(M,i);
    Ravg(:,i+M+1) = mean(real(R(:,idx)),2);
end
figure(111), imagesc(Ravg')
% r = [-abs_dist(end:-1:1) 0 abs_dist(1:end)];
% figure(63), plot(r,Ravg(round(Nlags/2),:))

end


