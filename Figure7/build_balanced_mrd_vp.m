function files = build_balanced_mrd_vp()
%BUILD_BALANCED_MRD_VP One balanced flow-MRD draw from the complete VP table.
%
% MRD-negative specimens contribute one virtual patient. MRD-positive
% specimens contribute four, and nine of them contribute a fifth, so the
% table has 85 rows in each class. Sampling is without replacement within
% a specimen.

thisDir = fileparts(mfilename('fullpath'));
outputDir = fullfile(thisDir, 'outputs');
if ~isfolder(outputDir), mkdir(outputDir); end
completeVPFile = fullfile(thisDir, 'Data', 'PB_BM_merged_with_MRDflow_BCR.csv');
assert(isfile(completeVPFile), 'Missing VP table: %s', completeVPFile);

featureNames = ["M11","M12","M13","M14","M21","M22","M23","M24", ...
    "M31","M32","M33","M34","M41","M42","M43","M44"];
T = readtable(completeVPFile, 'VariableNamingRule', 'preserve', 'TextType', 'string');
mrd = strtrim(string(T{:, 'MRD(flow)'}));
keep = ismember(mrd, ["Negative","Positive (< 10^-3)","Positive (> 10^-3)"]);
T = T(keep, :);
mrd = mrd(keep);
block = cumsum(double(T.VirtualPatient) == 1);
isPos = mrd == "Positive (> 10^-3)";
negBlocks = unique(block(~isPos), 'stable');
posBlocks = unique(block(isPos), 'stable');
assert(numel(negBlocks) == 85, 'Expected 85 MRD-negative specimens.');
assert(numel(posBlocks) == 19, 'Expected 19 MRD-positive specimens.');
for b = [negBlocks; posBlocks]'
    rows = find(block == b);
    assert(numel(rows) == 100, 'Specimen block %g does not have 100 rows.', b);
    assert(numel(unique(isPos(rows))) == 1, 'Specimen block %g has mixed MRD labels.', b);
end

rng(1729, 'twister');
selected = zeros(170, 1);
cursor = 0;
for b = negBlocks'
    rows = find(block == b);
    cursor = cursor + 1;
    selected(cursor) = rows(randi(numel(rows)));
end
extraBlocks = posBlocks(randperm(numel(posBlocks), 9));
for b = posBlocks'
    rows = find(block == b);
    k = 4 + ismember(b, extraBlocks);
    chosen = rows(randperm(numel(rows), k));
    selected(cursor + (1:k)) = chosen;
    cursor = cursor + k;
end
assert(cursor == 170);
selected = selected(randperm(numel(selected)));

Tout = T(selected, :);
sourceBlock = block(selected);
sourceRow = selected;
Tout = addvars(Tout, sourceBlock, sourceRow, 'Before', 1);
outFile = fullfile(outputDir, 'balanced_mrd_vp_one_draw.csv');
writetable(Tout, outFile);

manifest = table(sourceRow, sourceBlock, string(T.ID(selected)), ...
    double(T.VirtualPatient(selected)), isPos(selected), ...
    'VariableNames', {'sourceRow','sourceBlock','sourceID', ...
    'virtualPatient','mrdPositive'});
manifestFile = fullfile(outputDir, 'balanced_mrd_vp_draw_manifest.csv');
writetable(manifest, manifestFile);
fprintf('Balanced draw: %d MRD- and %d MRD+ rows, %d source specimens.\n', ...
    sum(~manifest.mrdPositive), sum(manifest.mrdPositive), ...
    numel(unique(manifest.sourceBlock)));
fprintf('Wrote %s\n', outFile);

files.table = string(outFile);
files.manifest = string(manifestFile);
end
