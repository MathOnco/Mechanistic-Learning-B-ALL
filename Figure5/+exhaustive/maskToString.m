function value = maskToString(mask, featureNames)
%MASKTOSTRING Convert one bit mask to comma-separated transition names.

cols = exhaustive.indicesFromMask(mask, numel(featureNames));
if isempty(cols), value = ""; else, value = strjoin(featureNames(cols), ','); end
end
