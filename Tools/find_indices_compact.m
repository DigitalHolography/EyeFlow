function linear_idx = find_indices_compact(M, s)
    % Compact version that returns only linear indices
    [i, j] = find(bsxfun(@minus, (1:M)', 1:M) == s);
    [sorted_i, idx] = sort(i);
    linear_idx = sub2ind([M, M], sorted_i, j(idx));
end