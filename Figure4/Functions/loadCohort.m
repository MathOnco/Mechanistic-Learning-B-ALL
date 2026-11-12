function C = loadCohort(cfg, outcome)
%LOADCOHORT Build one analysis cohort with all predictor views aligned.

Dm = readTransposed(cfg.markovFile);
Df = readTransposed(cfg.flowFile);
if Dm.nSamples ~= Df.nSamples
    error('mechlearn:SampleMismatch', 'markov.csv and flow.csv have different sample counts.');
end

markovNames = ["M11","M12","M13","M14","M21","M22","M23","M24", ...
    "M31","M32","M33","M34","M41","M42","M43","M44"];
flowNames = ["x1","x2","x3","x4"];

Xmarkov = nan(Dm.nSamples, numel(markovNames));
Xflow = nan(Df.nSamples, numel(flowNames));
for j = 1:numel(markovNames)
    Xmarkov(:, j) = getNumericRow(Dm, markovNames(j))';
end
for j = 1:numel(flowNames)
    Xflow(:, j) = getNumericRow(Df, flowNames(j))';
end

patientID = getStringRow(Dm, 'PatientID')';
sampleID = getStringRow(Dm, 'SampleID')';
age = getNumericRow(Dm, 'Age')';
sex = getStringRow(Dm, 'Sex')';
specimen = getStringRow(Dm, 'Specimen')';
progression = getStringRow(Dm, 'DiseaseProgression')';
bcr = getStringRow(Dm, 'BCR::ABL1/like')';
mrdFlow = getStringRow(Dm, 'MRD(flow)')';
mrdMolecular = getStringRow(Dm, 'MRD(molecular)')';
relapse = getStringRow(Dm, '3YrRelapse')';

clinicalNames = ["Age","Sex","Specimen","MDS_ICC","IKZF1", ...
    "CDKN2AB","RAS_RAF_MAPK_or_FLT3","TP53"];
clinicalRaw = table(age, sex, specimen, ...
    getStringRow(Dm, 'MDS/ICC')', ...
    getStringRow(Dm, 'IKZF1')', ...
    getStringRow(Dm, 'CDKN2A/B')', ...
    getStringRow(Dm, 'RAS/RAF/MAPK or FLT3')', ...
    getStringRow(Dm, 'TP53')', ...
    'VariableNames', cellstr(clinicalNames));

isDiagnosis = ~strcmpi(progression, 'Remission') & ~strcmpi(progression, 'Relapse');
switch lower(string(outcome))
    case "mrdflow"
        eligible = mrdFlow == "Positive (> 10^-3)" | mrdFlow == "Negative" | ...
            mrdFlow == "Positive (< 10^-3)";
        y = mrdFlow == "Positive (> 10^-3)";
        positiveLabel = "MRD+ (>10^-3)";
        negativeLabel = "MRD-";
    case "bcr"
        eligible = bcr == "BCR::ABL1" | bcr == "Negative";
        y = bcr == "BCR::ABL1";
        positiveLabel = "BCR::ABL1";
        negativeLabel = "Negative";
    case "mrdmolecular"
        eligible = mrdMolecular == "Positive (> 10^-3)" | mrdMolecular == "Negative" | ...
            mrdMolecular == "Positive (< 10^-3)";
        y = mrdMolecular == "Positive (> 10^-3)";
        positiveLabel = "Molecular MRD+ (>10^-3)";
        negativeLabel = "Molecular MRD-";
    case "relapse"
        eligible = relapse == "Yes" | relapse == "No";
        y = relapse == "Yes";
        positiveLabel = "3-year relapse";
        negativeLabel = "No 3-year relapse";
    otherwise
        error('mechlearn:UnknownOutcome', 'Unknown outcome: %s', outcome);
end

finitePredictors = all(isfinite(Xmarkov), 2) & all(isfinite(Xflow), 2);
keep = isDiagnosis & eligible & finitePredictors;

C.outcome = string(outcome);
C.y = double(y(keep));
C.patientID = patientID(keep);
C.sampleID = sampleID(keep);
C.specimen = specimen(keep);
C.progression = progression(keep);
C.bcr = bcr(keep);
C.mrdFlow = mrdFlow(keep);
C.mrdMolecular = mrdMolecular(keep);
C.relapse = relapse(keep);
C.Xmarkov = Xmarkov(keep, :);
C.markovNames = markovNames;
C.Xflow = Xflow(keep, :);
C.flowNames = flowNames;
C.clinicalRaw = clinicalRaw(keep, :);
C.positiveLabel = positiveLabel;
C.negativeLabel = negativeLabel;
C.originalIndex = find(keep);
end
