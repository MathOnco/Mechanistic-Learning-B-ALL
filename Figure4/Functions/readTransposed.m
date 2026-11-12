function D = readTransposed(filename)
%READTRANSPOSED Read a variables-by-rows CSV without naive comma splitting.

if exist(filename, 'file') ~= 2
    error('reviewer:MissingFile', 'Input file not found: %s', filename);
end
raw = readcell(filename, 'Delimiter', ',');
rowNames = strtrim(string(raw(:, 1)));
D.filename = filename;
D.raw = raw;
D.rowNames = rowNames;
D.nSamples = size(raw, 2) - 1;
end
