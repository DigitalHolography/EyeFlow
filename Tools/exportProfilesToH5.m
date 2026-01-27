function exportProfilesToH5(name, v_cell, v_safe_cell, v_profiles_cell)
% Function to save the velocities profiles 
    arguments
        name string
        v_cell cell
        v_safe_cell cell
        v_profiles_cell cell
    end

    ToolBox = getGlobalToolBox;
    params = ToolBox.getParams;
    sys_idx_list = ToolBox.Cache.sysIdxList;

    v_mat = toArray(v_cell);
    v_safe_mat = toArray(v_safe_cell);
    v_profiles_mat = toArray4D(v_profiles_cell);

    bandLimitedSignalHarmonicCount = params.json.PulseAnalysis.BandLimitedSignalHarmonicCount;

    %{
    [velocitySignalPerBeatPerSegment_whole, velocitySignalPerBeatPerSegmentFFT_whole, velocitySignalPerBeatPerSegmentBandLimited_whole] = perBeatSignalAnalysisMat(v_safe_mat, sys_idx_list, bandLimitedSignalHarmonicCount);
    velocitySignalPerBeatPerSegment_whole            = mat2cell4D_shape(velocitySignalPerBeatPerSegment_whole);
    velocitySignalPerBeatPerSegmentFFT_whole         = mat2cell4D_shape(velocitySignalPerBeatPerSegmentFFT_whole);
    velocitySignalPerBeatPerSegmentBandLimited_whole = mat2cell4D_shape(velocitySignalPerBeatPerSegmentBandLimited_whole);

    [velocitySignalPerBeatPerSegment_trunc, velocitySignalPerBeatPerSegmentFFT_trunc, velocitySignalPerBeatPerSegmentBandLimited_trunc] = perBeatSignalAnalysisMat(v_mat, sys_idx_list, bandLimitedSignalHarmonicCount);
    velocitySignalPerBeatPerSegment_trunc            = mat2cell4D_shape(velocitySignalPerBeatPerSegment_trunc);
    velocitySignalPerBeatPerSegmentFFT_trunc         = mat2cell4D_shape(velocitySignalPerBeatPerSegmentFFT_trunc);
    velocitySignalPerBeatPerSegmentBandLimited_trunc = mat2cell4D_shape(velocitySignalPerBeatPerSegmentBandLimited_trunc);
    %}

    % [profilePerBeatPerSegments, profilePerBeatPerSegmentsFFT, profilePerBeatPerSegmentsBandLimited] = perBeatProfileAnalysisMat(v_profiles_mat, sys_idx_list, bandLimitedSignalHarmonicCount);
    % profilePerBeatPerSegments            = mat2cell5D_shape(profilePerBeatPerSegments);
    % profilePerBeatPerSegmentsFFT         = mat2cell5D_shape(profilePerBeatPerSegmentsFFT);
    % profilePerBeatPerSegmentsBandLimited = mat2cell5D_shape(profilePerBeatPerSegmentsBandLimited);

    % ToolBox.Output.add("profilePerBeatPerSegments"            + capitalize(name), profilePerBeatPerSegments,            h5path = capitalize(name) + "/PerBeat/Segments/profilePerBeatPerSegments",            keepSize=true);
    % ToolBox.Output.add("profilePerBeatPerSegmentsFFT_abs"     + capitalize(name), abs(profilePerBeatPerSegmentsFFT),    h5path = capitalize(name) + "/PerBeat/Segments/profilePerBeatPerSegmentsFFT_abs",     keepSize=true);
    % ToolBox.Output.add("profilePerBeatPerSegmentsFFT_arg"     + capitalize(name), angle(profilePerBeatPerSegmentsFFT),  h5path = capitalize(name) + "/PerBeat/Segments/profilePerBeatPerSegmentsFFT_arg",     keepSize=true);
    % ToolBox.Output.add("profilePerBeatPerSegmentsBandLimited" + capitalize(name), profilePerBeatPerSegmentsBandLimited, h5path = capitalize(name) + "/PerBeat/Segments/profilePerBeatPerSegmentsBandLimited", keepSize=true);

    %{ 
    ToolBox.Output.add("velocitySignalPerBeatPerSegmentWhole"            + capitalize(name), velocitySignalPerBeatPerSegment_whole,            h5path = capitalize(name) + "/PerBeat/Segments/velocitySignalPerBeatPerSegmentWhole",            keepSize=true);
    ToolBox.Output.add("velocitySignalPerBeatPerSegmentWholeFFT_abs"     + capitalize(name), abs(velocitySignalPerBeatPerSegmentFFT_whole),    h5path = capitalize(name) + "/PerBeat/Segments/velocitySignalPerBeatPerSegmentWholeFFT_abs",     keepSize=true);
    ToolBox.Output.add("velocitySignalPerBeatPerSegmentWholeFFT_arg"     + capitalize(name), angle(velocitySignalPerBeatPerSegmentFFT_whole),  h5path = capitalize(name) + "/PerBeat/Segments/velocitySignalPerBeatPerSegmentWholeFFT_arg",     keepSize=true);
    ToolBox.Output.add("velocitySignalPerBeatPerSegmentWholeBandLimited" + capitalize(name), velocitySignalPerBeatPerSegmentBandLimited_whole, h5path = capitalize(name) + "/PerBeat/Segments/velocitySignalPerBeatPerSegmentWholeBandLimited", keepSize=true);

    ToolBox.Output.add("velocitySignalPerBeatPerSegmentTrunc"            + capitalize(name), velocitySignalPerBeatPerSegment_trunc,            h5path = capitalize(name) + "/PerBeat/Segments/velocitySignalPerBeatPerSegmentTrunc",            keepSize=true);
    ToolBox.Output.add("velocitySignalPerBeatPerSegmentTruncFFT_abs"     + capitalize(name), abs(velocitySignalPerBeatPerSegmentFFT_trunc),    h5path = capitalize(name) + "/PerBeat/Segments/velocitySignalPerBeatPerSegmentTruncFFT_abs",     keepSize=true);
    ToolBox.Output.add("velocitySignalPerBeatPerSegmentTruncFFT_arg"     + capitalize(name), angle(velocitySignalPerBeatPerSegmentFFT_trunc),  h5path = capitalize(name) + "/PerBeat/Segments/velocitySignalPerBeatPerSegmentTruncFFT_arg",     keepSize=true);
    ToolBox.Output.add("velocitySignalPerBeatPerSegmentTruncBandLimited" + capitalize(name), velocitySignalPerBeatPerSegmentBandLimited_trunc, h5path = capitalize(name) + "/PerBeat/Segments/velocitySignalPerBeatPerSegmentTruncBandLimited", keepSize=true);
    %}

    %perBeatSignalAnalysisMat(v_mat,         "Trunc", sys_idx_list, bandLimitedSignalHarmonicCount);
    ToolBox = getGlobalToolBox;

    [velocitySignalPerBeatPerSegment, velocitySignalPerBeatPerSegmentFFT, velocitySignalPerBeatPerSegmentBandLimited] = perBeatSignalAnalysisMat_handle(v_mat, sys_idx_list, bandLimitedSignalHarmonicCount);
    velocitySignalPerBeatPerSegment            = mat2cell4D_shape(velocitySignalPerBeatPerSegment);
    velocitySignalPerBeatPerSegmentFFT         = mat2cell4D_shape(velocitySignalPerBeatPerSegmentFFT);
    velocitySignalPerBeatPerSegmentBandLimited = mat2cell4D_shape(velocitySignalPerBeatPerSegmentBandLimited);

    ToolBox.Output.add("VelocitySignalPerBeatPerSegment"                 + capitalize(name), velocitySignalPerBeatPerSegment,            h5path = capitalize(name) + "/VelocityPerBeat/Segments/VelocitySignalPerBeatPerSegment",            keepSize=true);
    ToolBox.Output.add("VelocitySignalPerBeatPerSegment" + "FFT_abs"     + capitalize(name), abs(velocitySignalPerBeatPerSegmentFFT),    h5path = capitalize(name) + "/VelocityPerBeat/Segments/VelocitySignalPerBeatPerSegmentFFT_abs",     keepSize=true);
    ToolBox.Output.add("VelocitySignalPerBeatPerSegment" + "FFT_arg"     + capitalize(name), angle(velocitySignalPerBeatPerSegmentFFT),  h5path = capitalize(name) + "/VelocityPerBeat/Segments/VelocitySignalPerBeatPerSegmentFFT_arg",     keepSize=true);
    ToolBox.Output.add("VelocitySignalPerBeatPerSegment" + "BandLimited" + capitalize(name), velocitySignalPerBeatPerSegmentBandLimited, h5path = capitalize(name) + "/VelocityPerBeat/Segments/VelocitySignalPerBeatPerSegmentBandLimited", keepSize=true);

    out = mode1(velocitySignalPerBeatPerSegment);
    exportMode1Struct(out, name);

    % ToolBox.Output.add("velocity_trunc_seg_mean_" + name, v_mat,      h5path = capitalize(name) + "/CrossSections/velocity_trunc_seg_mean");
    ToolBox.Output.add("velocity_seg_mean_" + name, v_safe_mat, h5path = capitalize(name) + "/CrossSections/VelocityPerSegment", keepSize = true, unit = "mm/s", dimDesc = ["Branch", "Circle"]);

    ToolBox.Output.add("VelocityProfileSeg" + name, v_profiles_mat, h5path = capitalize(name) + "/CrossSections/VelocityProfileSeg", keepSize = true);
end

function [velocitySignalPerBeat, velocitySignalPerBeatFFT, velocitySignalPerBeatBandLimited] = perBeatSignalAnalysisMat_handle(v_mat, sys_idx_list, bandLimitedSignalHarmonicCount)
    arguments
        v_mat,
        sys_idx_list,
        bandLimitedSignalHarmonicCount
    end
    [c_size, b_size, ~] = size(v_mat);
    velocitySignalPerBeat            = cell(c_size, b_size);
    velocitySignalPerBeatFFT         = cell(c_size, b_size);
    velocitySignalPerBeatBandLimited = cell(c_size, b_size);

    for c_idx = 1:c_size
        for b_idx = 1:b_size
            [velocitySignalPerBeat{c_idx, b_idx}, velocitySignalPerBeatFFT{c_idx, b_idx}, velocitySignalPerBeatBandLimited{c_idx, b_idx}] = perBeatSignalAnalysis(v_mat(c_idx, b_idx, :), sys_idx_list, bandLimitedSignalHarmonicCount);
        end
    end
end

function [profilePerBeat, profilePerBeatFFT, profilePerBeatBandLimited] = perBeatProfileAnalysisMat(v_profiles_4d, sys_idx_list, bandLimitedSignalHarmonicCount)
    arguments
        v_profiles_4d
        sys_idx_list
        bandLimitedSignalHarmonicCount
    end

    [c_size, b_size, ~, ~] = size(v_profiles_4d);

    profilePerBeat            = cell(c_size, b_size);
    profilePerBeatFFT         = cell(c_size, b_size);
    profilePerBeatBandLimited = cell(c_size, b_size);

    for c_idx = 1:c_size
        for b_idx = 1:b_size

            profile = squeeze(v_profiles_4d(c_idx, b_idx, :, :));

            if all(isnan(profile), 'all')
                continue
            end

            [profilePerBeat{c_idx, b_idx}, profilePerBeatFFT{c_idx, b_idx}, profilePerBeatBandLimited{c_idx, b_idx}] = perBeatProfileAnalysis(profile, sys_idx_list, bandLimitedSignalHarmonicCount);
        end
    end
end

function v_array = toArray(v_cell)
    [rows, cols] = size(v_cell);

    firstIdx = find(~cellfun(@isempty, v_cell), 1);
    if isempty(firstIdx)
        error("The cell array is completely empty.");
    end

    vecLen = length(v_cell{firstIdx});
    
    v_array = nan(rows, cols, vecLen, "double");
    
    for r = 1:rows
        for c = 1:cols
            if ~isempty(v_cell{r,c})
                v_array(r, c, :) = double(v_cell{r,c});
            end
        end
    end
end

function v_4d = toArray4D(v_cell)
    [rows, cols] = size(v_cell);

    firstIdx = find(~cellfun(@isempty, v_cell), 1);
    if isempty(firstIdx)
        v_4d = []; return;
    end
    
    innerCell = v_cell{firstIdx};
    numNested = numel(innerCell);
    vecLen = numel(innerCell{1});

    v_4d = nan(rows, cols, numNested, vecLen, "double");

    for i = 1:numel(v_cell)
        if ~isempty(v_cell{i})
            [r, c] = ind2sub([rows, cols], i);
            
            nestedMatrix = vertcat(v_cell{i}{:});
            v_4d(r, c, :, :) = nestedMatrix;
        end
    end
end

function array = mat2cell4D_shape(input_cell)
    arguments
        input_cell cell 
    end

    [x, y] = size(input_cell);
    [z, a] = size(input_cell{1, 1});
    array = cat(4, input_cell{:});     % z × a × (x*y)
    array = reshape(array, z, a, x, y);
    array = permute(array, [3 4 1 2]); % x × y × z × a
end

function array5D = mat2cell5D_shape(input_cell)
    arguments
        input_cell cell
    end

    % grid size
    [C, B] = size(input_cell);

    % find first non-empty cell
    firstIdx = find(~cellfun(@isempty, input_cell), 1);
    if isempty(firstIdx)
        array5D = [];
        return
    end

    % infer dimensions from first valid cell
    sample = input_cell{firstIdx};   % [beats × Nfft × pixels]
    [numBeats, Nfft, numPixels] = size(sample);

    % preallocate output
    array5D = nan(C, B, numBeats, Nfft, numPixels, 'double');

    % fill array
    for c = 1:C
        for b = 1:B
            if ~isempty(input_cell{c,b})
                array5D(c, b, :, :, :) = input_cell{c,b};
            end
        end
    end
end

function out = mode1(v)
% mode-1 decomposition/reconstruction:
%   w(t,b,k,r) = mu(b,k,r) + a1(b,k,r) * u1(t)
% from SVD of DC-removed waveforms in v(t,b,k,r).
%
% INPUT
%   v : T x B x K x R  (may contain NaNs in mapping and/or time)
%
% OUTPUT
%   out.u1        : T x 1
%   out.mu        : B x K x R
%   out.a1        : B x K x R (NaN where invalid)
%   out.w         : T x B x K x R (NaN where invalid)
%   out.wc        : T x B x K x R (NaN where invalid)
%   out.validMask : B x K x R logical
%   out.s         : singular values

    v = permute(v, [4, 3, 2, 1]);

    if ndims(v) ~= 4
        error('Expected v to be 4D: (T x B x K x R). Got %s', mat2str(size(v)));
    end

    [T,B,K,R] = size(v);
    N = B*K*R;

    % --- 1) local mean over time (NaN-aware): mu is (1 x B x K x R) then squeeze -> (B x K x R)
    mu4 = mean(v, 1, 'omitnan');        % 1 x B x K x R
    mu  = squeeze(mu4);                % B x K x R

    % --- 2) center each waveform by its own mean
    vc = v - mu4;                      % implicit expansion (R2016b+)

    % --- 3) reshape to matrix Mc: (T x N)
    Mc_all = reshape(vc, T, N);

    % --- 4) keep only columns fully finite across time
    validCol = all(isfinite(Mc_all), 1);   % 1 x N
    Mc = Mc_all(:, validCol);              % T x Nvalid
    if size(Mc,2) < 2
        error('Not enough valid waveforms for SVD (Nvalid=%d).', size(Mc,2));
    end

    % --- 5) SVD and mode-1 scores
    [U,S,V] = svd(Mc, 'econ');
    s = diag(S);
    u1 = U(:,1);                 % T x 1
    a1_valid = s(1) * V(:,1);    % Nvalid x 1

    % --- 6) scatter scores back into full (1 x N), then reshape to (B x K x R)
    a1_all = nan(1, N);
    a1_all(validCol) = a1_valid;
    a1 = reshape(a1_all, [B, K, R]);

    validMask = reshape(validCol, [B, K, R]);

    % --- 7) reconstruct: wc = u1 * a1, then add mu
    % Build wc as (T x N), then reshape to (T x B x K x R)
    wc_all = u1 * a1_all;                 % (T x 1) * (1 x N) -> T x N
    wc = reshape(wc_all, [T, B, K, R]);   % T x B x K x R

    % Add back mu: need mu4 (1 x B x K x R) broadcast along time
    w = wc + mu4;                         % T x B x K x R

    % If you want invalid entries to remain NaN (they will, because a1_all is NaN there)
    % but ensure mu isn't accidentally finite there:
    % (Optional) force invalid to NaN:
    % w(:, ~validMask)  = NaN;  % not allowed directly for 4D; keep as-is unless needed

    mu = permute(mu, [3, 2, 1]);
    a1 = permute(a1, [3, 2, 1]);
    wc = permute(wc, [4, 3, 2, 1]);
    w = permute(w, [4, 3, 2, 1]);
    validMask = permute(validMask, [3, 2, 1]);

    out = struct();
    out.u1        = u1;
    out.mu        = mu;
    out.a1        = a1;
    out.wc        = wc;
    out.w         = w;
    out.validMask = validMask;
    out.s         = s;
end

function exportMode1Struct(out, name)
    arguments
        out
        name
    end

    ToolBox = getGlobalToolBox;

    var_name = "VelocityMode1SignalPerBeatPerSegment";
    h5path   = capitalize(name) + "/VelocityPerBeat/Segments/VelocityMode1SignalPerBeatPerSegment";

    ToolBox.Output.add(var_name + "_mu_"        + capitalize(name), out.mu, ...
        h5path      = h5path + "/mu",   ...
        dimDesc     = ["Beat", "Branch", "Circle"], ...
        keepSize    = true);
    ToolBox.Output.add(var_name + "_a1_"        + capitalize(name), out.a1, ...
        h5path      = h5path + "/a1", ...
        dimDesc     = ["Beat", "Branch", "Circle"], ...
        keepSize    = true);
    ToolBox.Output.add(var_name + "_u1_"        + capitalize(name), out.u1, ...
        h5path      = h5path + "/u1", ...
        dimDesc     = ["Time"], ...
        keepSize    = true);
    ToolBox.Output.add(var_name + "_wc_"        + capitalize(name), out.wc, ...
        h5path      = h5path + "/wc", ...
        dimDesc     = ["Time", "Beat", "Branch", "Circle"], ...
        keepSize    = true);
    ToolBox.Output.add(var_name + "_w_"         + capitalize(name), out.w, ...
        h5path      = h5path + "/w", ...
        dimDesc     = ["Time", "Beat", "Branch", "Circle"], ...
        keepSize    = true);
    ToolBox.Output.add(var_name + "_validMask_" + capitalize(name), out.validMask, ...
        h5path      = h5path + "/validMask", ...
        dimDesc     = ["Beat", "Branch", "Circle"], ...
        keepSize    = true);
    ToolBox.Output.add(var_name + "_s_"         + capitalize(name), out.s, ...
        h5path      = h5path + "/s", ...
        dimDesc     = ["Time"], ...
        keepSize    = true);
end

% +============================================+ %
% |                   OLD CODE                 | %
% +============================================+ %

%{
function exportProfilesToH5(name, maskLabel, v_cell, dv_cell)

% ToolBox = getGlobalToolBox;

% if nargin < 4
%     dv_cell = [];
% end

% [numCircles, numBranches] = size(maskLabel);
% % [numCircles, numBranches] = size(v_profiles_cell);

% numFrames = 0;
% i = 1;

% while numFrames <= 0
%     numFrames = size(v_cell{i}, 2);
%     i = i + 1;

%     if i > size(v_cell, 1) * size(v_cell, 2)
%         warning("Velocity profiles cells are all empty.")
%         break
%     end

% end

% L = zeros(size(maskLabel{1, 1}), 'uint8');
% idx = 1;
% % Process each circle and branch

% ProfilesTimeBranches = zeros([numFrames, numCircles, numBranches]);

% for cIdx = 1:numCircles

%     for bIdx = 1:numBranches

%         if ~isempty(v_cell{cIdx, bIdx})
%             current_v = v_cell{cIdx, bIdx};
%             cv = zeros(size(current_v{1}, 2), numFrames, "single");

%             for ff = 1:numFrames
%                 cv(:, ff) = current_v{ff};
%             end

%             if ~isempty(dv_cell)
%                 current_dv = dv_cell{cIdx, bIdx};
%                 cdv = zeros(size(current_dv{1}, 2), numFrames, "single");

%                 for ff = 1:numFrames
%                     cdv(:, ff) = current_dv{ff};
%                 end

%             else
%                 cdv = [];
%             end

%             L(maskLabel{cIdx, bIdx}) = idx;
%             ProfilesTimeBranches(:, cIdx, bIdx) = mean(cv, 1, 'omitnan');
%             ToolBox.Output.Extra.add(sprintf("Segments/%s_idx%d_c%d_b%d_VelocityProfiles", name, idx, cIdx, bIdx), cv, cdv, "mm/s");
%             idx = idx + 1;
%         end

%     end

% end

% ToolBox.Output.Extra.add(sprintf("Segments/%s_Segments_Labels", name), L);

end
%}
