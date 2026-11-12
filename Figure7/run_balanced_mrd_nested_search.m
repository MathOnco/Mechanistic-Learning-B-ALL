function run_balanced_mrd_nested_search()
%RUN_BALANCED_MRD_NESTED_SEARCH Nested best-k search on one balanced MRD draw.
%
% Penalties are the mean lambdas from the full virtual-patient nested CV.
% Subset choice and outer-test accuracy use 10 outer and 10 inner folds.
% All rows from one source specimen stay in the same fold.

thisDir = fileparts(mfilename('fullpath'));
addpath(thisDir);
addpath(fullfile(thisDir, 'Functions'));
outputDir = fullfile(thisDir, 'outputs');
if ~isfolder(outputDir), mkdir(outputDir); end

build_balanced_mrd_vp();
lambdaFile = fullfile(outputDir, 'vp_mrd_mean_lambdas.csv');
assert(isfile(lambdaFile), ...
    'Missing %s. Run collect_vp_mrd_mean_lambdas first.', lambdaFile);
lambdas = readtable(lambdaFile, 'TextType', 'string');

cfg = exhaustive_config();
cfg.outputDir = outputDir;
cfg.outcomes = "MRDflow";
cfg.models = ["ridge","lasso","elasticnet"];
cfg.views = "markov";
cfg.nRepeats = 1;
cfg.outerFolds = 10;
cfg.innerFolds = 10;
cfg.maxK = 16;
cfg.useParallel = true;
cfg.resumeFromCheckpoint = true;
warning('off', 'stats:lasso:MaxIterReached');
warning('off', 'stats:glmfit:IterationLimit');

if license('test', 'Distrib_Computing_Toolbox')
    functionsDir = fileparts(which('runEndpointModel'));
    pool = gcp('nocreate');
    if isempty(pool), pool = parpool; end
    pctRunOnAll(['addpath(''' strrep(functionsDir, '''', '''''') ''')']);
    fprintf('Parallel pool: %d workers\n', pool.NumWorkers);
else
    cfg.useParallel = false;
end

C = loadBalancedCohort(fullfile(outputDir, 'balanced_mrd_vp_one_draw.csv'));
fprintf(['Balanced MRD VP: n=%d (+%d/-%d), %d specimens. ' ...
    'Folds are assigned by specimen.\n'], ...
    numel(C.y), sum(C.y == 1), sum(C.y == 0), numel(unique(C.group)));

for m = 1:numel(cfg.models)
    model = cfg.models(m);
    hit = lambdas.model == model;
    assert(sum(hit) == 1, 'Expected one mean lambda for %s.', model);
    lambda = lambdas.meanLambda(hit);
    fprintf('\n%s fixed mean lambda %.15g\n', model, lambda);
    R = runEndpointModel(C, "MRDflow", model, lambda, cfg);
    saveResults(R, C, cfg);
end
fprintf('Balanced MRD nested search complete: %s\n', outputDir);
end

function C = loadBalancedCohort(filename)
featureNames = ["M11","M12","M13","M14","M21","M22","M23","M24", ...
    "M31","M32","M33","M34","M41","M42","M43","M44"];
T = readtable(filename, 'VariableNamingRule', 'preserve', 'TextType', 'string');
mrd = strtrim(string(T{:, 'MRD(flow)'}));
C.outcome = "MRDflow";
C.view = "markov";
C.X = double(T{:, cellstr(featureNames)});
C.y = double(mrd == "Positive (> 10^-3)");
C.group = double(T.sourceBlock);
C.patientID = string(T.ID);
C.sampleID = string(T.ID) + "__VP" + string(T.VirtualPatient);
C.featureNames = featureNames;
C.positiveLabel = "MRD+ (>10^-3)";
C.negativeLabel = "MRD-";
C.originalIndex = double(T.sourceRow);
assert(all(isfinite(C.X), 'all'), 'Balanced VP matrix has nonfinite values.');
assert(sum(C.y == 1) == 85 && sum(C.y == 0) == 85, ...
    'The balanced draw must contain 85 rows in each class.');
end
