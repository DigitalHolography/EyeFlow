function s = MEdisp(ME, path)

s = sprintf("==================================\nERROR\n==================================\n");
s = s + sprintf('Error while processing : %s\n', path);
s = s + sprintf("%s\n", ME.identifier);
s = s + sprintf("%s\n", ME.message);

for stackIdx = 1:size(ME.stack, 1)
    s = s + sprintf("%s : %s, line : %d\n", ME.stack(stackIdx).file, ME.stack(stackIdx).name, ME.stack(stackIdx).line);
end

if ME.identifier == "MATLAB:audiovideo:VideoReader:FileNotFound"

    s = s + sprintf("\nNo Raw File was found\n");

end

s = s + sprintf("==================================\n");

if nargout == 0
    fprintf(2,"%s", s);
end

end
