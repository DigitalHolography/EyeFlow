function [passed, error_report] = checkStructConsistency(in_cells)
% CHECKSTRUCTCONSISTENCY Checks if all non-empty structures in a cell array
% have identical field structures (including recursively).
%
% Usage:
%   [ok, msg] = checkStructConsistency(my_cells)
%
% Outputs:
%   passed       : Boolean (true if consistent, false otherwise)
%   error_report : String describing the first mismatch found

    passed = true;
    error_report = 'All structures are consistent.';

    % Flatten for easier looping
    in_cells = in_cells(:);
    
    % 1. Find the Reference (First non-empty cell)
    mask_not_empty = ~cellfun(@isempty, in_cells);
    if ~any(mask_not_empty)
        error_report = 'Cell array is completely empty.';
        return; 
    end
    
    idx_filled = find(mask_not_empty);
    ref_idx = idx_filled(1);
    ref_struct = in_cells{ref_idx};
    
    % Ensure the reference is actually a struct
    if ~isstruct(ref_struct)
        passed = false;
        error_report = sprintf('Reference element at index %d is not a struct (it is %s).', ...
            ref_idx, class(ref_struct));
        return;
    end

    % 2. Compare every other non-empty cell against the reference
    for i = 2:length(idx_filled)
        curr_idx = idx_filled(i);
        curr_struct = in_cells{curr_idx};
        
        % Start recursive comparison
        [is_match, fail_reason] = compare_recursive(ref_struct, curr_struct, 'root');
        
        if ~is_match
            passed = false;
            error_report = sprintf('Mismatch at Index %d: %s', curr_idx, fail_reason);
            return; % Stop at first error
        end
    end
end

function [ok, msg] = compare_recursive(ref, target, path)
    % Compare two variables. If they are structs, check fields and recurse.
    ok = true;
    msg = '';

    % 1. Check if both are structs
    if isstruct(ref) && ~isstruct(target)
        ok = false;
        msg = sprintf('At "%s": Expected struct, got %s.', path, class(target));
        return;
    elseif ~isstruct(ref) && isstruct(target)
        ok = false;
        msg = sprintf('At "%s": Expected %s, got struct.', path, class(ref));
        return;
    elseif ~isstruct(ref) && ~isstruct(target)
        % Neither are structs, we don't care about values, just structure. 
        % So we are done here.
        return;
    end

    % If we are here, BOTH are structs.
    % Note: If ref is a struct ARRAY (1xN), all elements share the same fields.
    % We only need to check the first element of the array for field names.
    if numel(ref) > 0
        ref = ref(1); 
    end
    if numel(target) > 0
        target = target(1);
    end

    % 2. Compare Field Names
    fields_ref = sort(fieldnames(ref));
    fields_tgt = sort(fieldnames(target));

    if ~isequal(fields_ref, fields_tgt)
        ok = false;
        
        % Find exactly what is wrong for the message
        missing = setdiff(fields_ref, fields_tgt);
        extra = setdiff(fields_tgt, fields_ref);
        
        detail = '';
        if ~isempty(missing)
            detail = [detail ' Missing: [' strjoin(missing', ', ') ']'];
        end
        if ~isempty(extra)
            detail = [detail ' Extra: [' strjoin(extra', ', ') ']'];
        end
        
        msg = sprintf('Fields differ at "%s".%s', path, detail);
        return;
    end

    % 3. Recursively check contents of fields
    for k = 1:numel(fields_ref)
        fn = fields_ref{k};
        new_path = [path '.' fn];
        
        [sub_ok, sub_msg] = compare_recursive(ref.(fn), target.(fn), new_path);
        
        if ~sub_ok
            ok = false;
            msg = sub_msg;
            return;
        end
    end
end