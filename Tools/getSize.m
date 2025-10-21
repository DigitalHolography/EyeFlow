function getSize(this)
%GETSIZE Get the size of the object in bytes
props = properties(this);
totSize = 0;

for ii = 1:length(props)
    currentProperty = getfield(this, char(props(ii)));
    s = whos('currentProperty');
    totSize = totSize + s.bytes;
end

if totSize < 1024
    fprintf(1, '%d bytes\n', totSize);
elseif totSize < 1024 ^ 2
    totSize = totSize / 1024;
    fprintf(1, '%.2f KB\n', totSize);
elseif totSize < 1024 ^ 3
    totSize = totSize / 1024 ^ 2;
    fprintf(1, '%.2f MB\n', totSize);
else
    totSize = totSize / 1024 ^ 3;
    fprintf(1, '%.2f GB\n', totSize);
end

end
