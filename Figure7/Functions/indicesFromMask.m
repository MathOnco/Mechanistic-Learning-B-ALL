function cols = indicesFromMask(mask, p)
%INDICESFROMMASK Decode one uint16 feature mask.

cols = find(bitget(uint16(mask), 1:p) ~= 0);
end
