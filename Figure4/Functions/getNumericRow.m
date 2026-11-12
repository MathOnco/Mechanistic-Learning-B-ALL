function row = getNumericRow(D, name)
%GETNUMERICROW Retrieve a numeric row from a transposed data file.

s = getStringRow(D, name);
row = str2double(s);
end
