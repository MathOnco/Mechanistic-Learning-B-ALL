function cfg = _config()


thisDir = fileparts(mfilename('fullpath'));
cfg.analysisDir = thisDir;
cfg.dataDir = fullfile(thisDir, 'Data');

cfg.smokeMode = strcmp(getenv('REVIEWER_SMOKE_MODE'), '1');
if cfg.smokeMode
    cfg.outputDir = fullfile(thisDir, 'smoke_test_outputs');
else
    cfg.outputDir = fullfile(thisDir, 'outputs');
end
if ~exist(cfg.outputDir, 'dir')
    mkdir(cfg.outputDir);
end

cfg.markovFile = fullfile(cfg.dataDir, 'markov.csv');
cfg.flowFile = fullfile(cfg.dataDir, 'flow.csv');
cfg.virtualPatientFile = fullfile(cfg.dataDir, 'merged_markov_training.csv');

cfg.nRepeats = 100;
cfg.outerFolds = 5;
cfg.innerFolds = 5;
cfg.bootstrapIterations = 2000;
cfg.randomSeed = 1729;

% Hyperparameters selected only on each outer-training partition.
cfg.lambdaGrid = logspace(-4, 2, 13);
cfg.elasticNetAlpha = 0.5;
cfg.rfMinLeafGrid = [1 3 5 10];
cfg.numTrees = 300;

% Primary metrics are AUC and balanced accuracy. Accuracy is retained only
% for comparability with the submitted manuscript; no composite metric is
% used for tuning or model selection.
cfg.primaryTuneMetric = "auc";
cfg.models = ["ridge", "lasso", "elasticnet", "randomforest"];
cfg.views = ["markov", "flow", "clinical_genomic", "markov_plus_clinical"];

% Source-state abundance thresholds used in the sensitivity analysis.
cfg.abundanceThresholds = [0 0.005 0.01 0.025 0.05];

% VP augmentation-size analysis. These rows intentionally remain the unit
% of analysis, matching the authors' stated modeling objective.
cfg.vpCounts = [1 5 10 25 50 100];
cfg.vpRegimeRepeats = 30;

if cfg.smokeMode
    cfg.nRepeats = 2;
    cfg.bootstrapIterations = 100;
    cfg.lambdaGrid = [0.1 1 10];
    cfg.rfMinLeafGrid = [1 5];
    cfg.numTrees = 30;
    cfg.models = ["ridge", "lasso"];
    cfg.views = ["markov", "flow"];
    cfg.vpRegimeRepeats = 2;
end
end
