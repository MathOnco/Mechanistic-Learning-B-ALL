function D = readTransposed(filename)
%READTRANSPOSED Read a variables-by-rows CSV without naive comma splitting.

if exist(filename, 'file') ~= 2
    error('exhaustive:MissingFile', 'Input file not found: %s', filename);
end
raw = readcell(filename, 'Delimiter', ',');
D.filename = filename;
D.raw = raw;
D.rowNames = strtrim(string(raw(:, 1)));
D.nSamples = size(raw, 2) - 1;
end
