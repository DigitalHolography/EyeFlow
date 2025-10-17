function test_disk_speed(path)
    % TEST_DISK_SPEED  Measure write and read speeds at a given path.
    %
    % Usage:
    %   test_disk_speed('C:\temp')
    %   test_disk_speed('/mnt/data')
    %
    % Example output:
    %   Writing 500 MB...
    %   Write speed: 320.45 MB/s
    %   Reading 500 MB...
    %   Read speed: 355.12 MB/s

    if nargin < 1
        error('Please specify a path, e.g. test_disk_speed(''C:\temp'')');
    end

    % Parameters
    fileSizeMB = 500;                % File size in MB (adjust as needed)
    fileName = fullfile(path, 'disk_speed_test.bin');
    blockSize = 1024 * 1024;         % 1 MB blocks
    numBlocks = fileSizeMB;

    % Generate random data block
    dataBlock = rand(1, blockSize, 'single'); 

    % --- WRITE TEST ---
    fprintf('Writing %d MB to %s...\n', fileSizeMB, path);
    fid = fopen(fileName, 'w');
    if fid == -1
        error('Cannot open file for writing: %s', fileName);
    end
    tic;
    for i = 1:numBlocks
        fwrite(fid, dataBlock, 'single');
    end
    fclose(fid);
    writeTime = toc;
    writeSpeed = fileSizeMB / writeTime;
    fprintf('Write speed: %.2f MB/s\n', writeSpeed);

    % --- READ TEST ---
    fprintf('Reading %d MB...\n', fileSizeMB);
    fid = fopen(fileName, 'r');
    if fid == -1
        error('Cannot open file for reading: %s', fileName);
    end
    tic;
    for i = 1:numBlocks
        fread(fid, [1, blockSize], 'single');
    end
    fclose(fid);
    readTime = toc;
    readSpeed = fileSizeMB / readTime;
    fprintf('Read speed: %.2f MB/s\n', readSpeed);

    % --- CLEANUP ---
    delete(fileName);
    fprintf('Temporary file deleted.\n');
end
