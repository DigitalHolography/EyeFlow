function warning_s(varargin)
    % warning_s Issues a warning with only the top-level stack trace.
    % Usage mimics standard warning():
    %   warningTopStack('Message')
    %   warningTopStack('Message', arg1, arg2)
    %   warningTopStack('ID:IdString', 'Message', ...)

    if nargin == 0
        return;
    end

    st = dbstack;
    if length(st) > 1
        caller = st(2);

        fullPath = which(caller.file);
        if isempty(fullPath)
            fullPath = caller.file;
        end
        
        % Escape single quotes in the path 
        % If path is "C:\User's\Code", we need "C:\User''s\Code" for the eval string
        safePath = strrep(fullPath, '''', '''''');

        linkStr = sprintf('<a href="matlab:opentoline(''%s'',%d)">%s (line %d)</a>', ...
            safePath, caller.line, caller.name, caller.line);
        locMsg = sprintf('\n> In %s', linkStr);
    else
        locMsg = '';
    end

    oldState = warning('off', 'backtrace');
    restoreState = onCleanup(@() warning(oldState));

    arg1 = varargin{1};

    isText = ischar(arg1) || (isstring(arg1) && isscalar(arg1));
    
    hasID = false;
    if isText && nargin > 1
        % Heuristic: IDs contain ':' and usually no spaces (e.g., 'Msg:Id')
        txt = char(arg1); 
        if contains(txt, ':') && ~contains(txt, ' ')
            hasID = true;
        end
    end

    if hasID
        % Structure: warning(ID, Message, Args...)
        id = varargin{1};       % Keep original type (string or char)
        rawMsg = varargin{2};   % Message template
        args = varargin(3:end); % Formatting arguments
        
        warning(id, [char(rawMsg) '%s'], args{:}, locMsg);
    else
        % Structure: warning(Message, Args...)
        rawMsg = varargin{1};
        args = varargin(2:end);
        
        finalMsg = [char(rawMsg), locMsg];
        
        warning([char(rawMsg) '%s'], args{:}, locMsg);
    end
end