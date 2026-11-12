function [X, names] = encodeView(C, viewName)
%ENCODEVIEW Construct a numeric predictor matrix for a named data view.

viewName = lower(string(viewName));
[Xclinical, clinicalNames] = encodeClinical(C.clinicalRaw);

switch viewName
    case "markov"
        X = C.Xmarkov;
        names = C.markovNames;
    case "flow"
        X = C.Xflow;
        names = C.flowNames;
    case "clinical_genomic"
        X = Xclinical;
        names = clinicalNames;
    case "markov_plus_clinical"
        X = [C.Xmarkov, Xclinical];
        names = [C.markovNames, clinicalNames];
    case "flow_plus_clinical"
        X = [C.Xflow, Xclinical];
        names = [C.flowNames, clinicalNames];
    case "key_markov"
        if lower(C.outcome) == "bcr"
            idx = C.markovNames == "M31";
        else
            idx = ismember(C.markovNames, ["M11","M14"]);
        end
        X = C.Xmarkov(:, idx);
        names = C.markovNames(idx);
    otherwise
        error('reviewer:UnknownView', 'Unknown predictor view: %s', viewName);
end
end

function [X, names] = encodeClinical(T)
% Excludes WHO Category and Risk Group because they can encode the target.
X = T.Age;
names = "Age";
categoricalVars = setdiff(string(T.Properties.VariableNames), "Age", 'stable');
for j = 1:numel(categoricalVars)
    varName = categoricalVars(j);
    s = strtrim(string(T.(varName)));
    s(ismissing(s) | s == "" | strcmpi(s, "NaN")) = "Missing";
    levels = unique(s, 'stable');
    % Reference coding avoids an exact dummy-variable dependency.
    for k = 2:numel(levels)
        X(:, end+1) = double(s == levels(k)); %#ok<AGROW>
        safeLevel = matlab.lang.makeValidName(char(levels(k)));
        names(end+1) = varName + "_" + string(safeLevel); %#ok<AGROW>
    end
end
end
