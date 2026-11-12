function row = getNumericRow(D, name)
%GETNUMERICROW Retrieve one numeric row from a transposed data file.

row = str2double(getStringRow(D, name));
end
