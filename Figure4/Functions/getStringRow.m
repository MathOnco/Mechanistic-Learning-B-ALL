function row = getStringRow(D, name)
%GETSTRINGROW Retrieve and normalize a row from a transposed data file.

idx = find(D.rowNames == string(name), 1);
if isempty(idx)
    error('reviewer:MissingVariable', 'Variable "%s" not found in %s.', name, D.filename);
end
cells = D.raw(idx, 2:end);
row = strings(1, numel(cells));
for i = 1:numel(cells)
    value = cells{i};
    if isempty(value)
        row(i) = "NaN";
    elseif ischar(value) || isstring(value)
        converted = strtrim(string(value));
        if isempty(converted) || all(ismissing(converted) | converted == "")
            row(i) = "NaN";
        else
            row(i) = converted(1);
        end
    elseif isnumeric(value) || islogical(value)
        if isscalar(value)
            row(i) = string(value);
        else
            row(i) = "NaN";
        end
    else
        converted = strtrim(string(value));
        if isempty(converted) || all(ismissing(converted) | converted == "")
            row(i) = "NaN";
        else
            row(i) = converted(1);
        end
    end
end
row(ismissing(row) | row == "") = "NaN";
end
