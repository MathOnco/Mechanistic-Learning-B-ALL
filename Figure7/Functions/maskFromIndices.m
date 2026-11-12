function mask = maskFromIndices(cols)
%MASKFROMINDICES Encode transition indices 1..16 in one uint16 bit mask.

mask = uint16(0);
for j = 1:numel(cols)
    mask = bitset(mask, cols(j), 1);
end
end
