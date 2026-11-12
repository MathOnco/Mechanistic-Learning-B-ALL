function C = loadCohort(cfg, outcome, view)
%LOADCOHORT Load one Markov-transition or Flow-abundance diagnosis cohort.

view = lower(string(view));
Dm = readTransposed(cfg.markovFile);
switch view
    case "markov"
        Dx = Dm;
        featureNames = cfg.markovFeatureNames;
    case "flow"
        Dx = readTransposed(cfg.flowFile);
        featureNames = cfg.flowFeatureNames;
        if Dm.nSamples ~= Dx.nSamples
            error('exhaustive:SampleMismatch', ...
                'markov.csv and flow.csv have different sample counts.');
        end
    otherwise
        error('exhaustive:UnknownView', 'Unknown predictor view: %s', view);
end

X = nan(Dx.nSamples, numel(featureNames));
for j = 1:numel(featureNames)
    X(:, j) = getNumericRow(Dx, featureNames(j))';
end

patientID = getStringRow(Dm, 'PatientID')';
sampleID = getStringRow(Dm, 'SampleID')';
progression = getStringRow(Dm, 'DiseaseProgression')';
bcr = getStringRow(Dm, 'BCR::ABL1/like')';
mrdFlow = getStringRow(Dm, 'MRD(flow)')';

isDiagnosis = ~strcmpi(progression, 'Remission') & ...
    ~strcmpi(progression, 'Relapse');
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
    otherwise
        error('exhaustive:UnknownOutcome', 'Unknown outcome: %s', outcome);
end

keep = isDiagnosis & eligible & all(isfinite(X), 2);
C.outcome = string(outcome);
C.view = view;
C.X = double(X(keep, :));
C.y = double(y(keep));
C.patientID = patientID(keep);
C.sampleID = sampleID(keep);
C.featureNames = featureNames;
C.positiveLabel = positiveLabel;
C.negativeLabel = negativeLabel;
C.originalIndex = find(keep);
end
