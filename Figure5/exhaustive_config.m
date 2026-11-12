function cfg = exhaustive_config()
%EXHAUSTIVE_CONFIG Configuration for the separate exhaustive nested analysis.

thisDir = fileparts(mfilename('fullpath'));
cfg.analysisDir = thisDir;
cfg.dataDir = fullfile(thisDir, 'Data');
cfg.markovFile = fullfile(cfg.dataDir, 'markov.csv');
cfg.flowFile = fullfile(cfg.dataDir, 'flow.csv');

cfg.smokeMode = strcmp(getenv('EXHAUSTIVE_SMOKE_MODE'), '1');
if cfg.smokeMode
    cfg.outputDir = fullfile(thisDir, 'smoke_test_outputs');
else
    cfg.outputDir = fullfile(thisDir, 'outputs_10fold');
end
if ~isfolder(cfg.outputDir), mkdir(cfg.outputDir); end

% One ten-fold outer evaluation. The inner ten-fold partition performs
% subset selection without using any outer-test observations.
cfg.nRepeats = 1;
cfg.outerFolds = 10;
cfg.innerFolds = 10;
cfg.randomSeed = 1729;

cfg.outcomes = ["MRDflow", "BCR"];
cfg.models = ["ridge", "lasso", "elasticnet"];
cfg.views = "markov";
cfg.elasticNetAlpha = 0.5;

% Equal-weight discrimination score used for exhaustive selection.
% Ordinary accuracy remains an output but has zero weight in selection.
% Other supported selectionMetric values are "auc", "balancedaccuracy",
% and "accuracy".
cfg.selectionMetric = "combo";
cfg.comboWeights = [0.00, 0.50, 0.50]; % accuracy, balanced accuracy, AUC
cfg.tieTolerance = 1e-12;

cfg.markovFeatureNames = ["M11","M12","M13","M14", ...
    "M21","M22","M23","M24", ...
    "M31","M32","M33","M34", ...
    "M41","M42","M43","M44"];
cfg.flowFeatureNames = ["x1","x2","x3","x4"];
cfg.maxK = Inf;

% Mean lambdas across the 500 outer-fold choices in the classifier benchmark.
cfg.lambdaTable = table( ...
    ["markov"; "markov"; "markov"; "markov"; "markov"; "markov"; ...
     "flow"; "flow"; "flow"; "flow"; "flow"; "flow"], ...
    ["MRDflow"; "MRDflow"; "MRDflow"; "BCR"; "BCR"; "BCR"; ...
     "MRDflow"; "MRDflow"; "MRDflow"; "BCR"; "BCR"; "BCR"], ...
    ["ridge"; "lasso"; "elasticnet"; "ridge"; "lasso"; "elasticnet"; ...
     "ridge"; "lasso"; "elasticnet"; "ridge"; "lasso"; "elasticnet"], ...
    [15.6948594699224; 0.025129399895361; 0.0489558106388473; ...
     45.0032666062912; 0.0817761452253173; 0.196806670559214; ...
     8.70313263894129; 0.0109903268158237; 0.0194097088637811; ...
     63.5554207865417; 0.0758708604466219; 0.206238186309145], ...
    'VariableNames', {'view','outcome','model','lambda'});

% EBIC is descriptive here; inner-CV performance selects the best subset.
cfg.ebicGamma = [0.5, 1.0];

cfg.useParallel = true;
cfg.resumeFromCheckpoint = true;
cfg.saveCheckpointEveryK = true;

if cfg.smokeMode
    cfg.outcomes = "MRDflow";
    cfg.models = "ridge";
    cfg.views = "markov";
    cfg.outerFolds = 2;
    cfg.innerFolds = 2;
    cfg.maxK = 2;
    cfg.useParallel = false;
end
end
