function setupParpool(requestedWorkers)
%SETUPPARPOOL Set up parallel pool with specified number of workers
%   SETUPPARPOOL(requestedWorkers) checks if a parallel pool exists and
%   creates one with the specified number of workers if needed.
%
%   INPUT:
%   requestedWorkers - Number of workers desired for parallel processing
%
%   EXAMPLE:
%   setupParpool(4) % Sets up a parallel pool with 4 workers

fprintf("\n----------------------------------\n");
fprintf('---- Setting up parallel pool ----\n');
fprintf("----------------------------------\n");

% Validate input
if nargin < 1
    error('Number of workers must be specified');
end

if ~isnumeric(requestedWorkers) || requestedWorkers < 1 || mod(requestedWorkers, 1) ~= 0
    error('Number of workers must be a positive integer');
end

% Check if parallel computing toolbox is available
if ~license('test', 'Distrib_Computing_Toolbox')
    warning('Parallel Computing Toolbox is not available');
    return;
end

% Check for existing parallel pool
poolobj = gcp('nocreate');

if isempty(poolobj)

    parpool(requestedWorkers); % create a new pool
    fprintf('Created new parallel pool with %d workers\n', requestedWorkers);

elseif poolobj.NumWorkers ~= requestedWorkers

    delete(poolobj); %close the current pool to create a new one with correct num of workers
    parpool(requestedWorkers);
    fprintf('Recreated parallel pool with %d workers (was %d)\n', ...
        requestedWorkers, poolobj.NumWorkers);

else

    fprintf('Parallel pool already exists with %d workers\n', requestedWorkers);

end

end
