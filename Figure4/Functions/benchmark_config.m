function cfg = benchmark_config()
%BENCHMARK_CONFIG Settings for the self-contained classifier benchmark.
%
% Data and outputs stay inside Final Files/Figure4. Set the environment

thisDir = fileparts(mfilename('fullpath'));
figureDir = fileparts(thisDir);
cfg.analysisDir = figureDir;
cfg.dataDir = fullfile(figureDir, 'Data');
cfg.smokeMode = strcmp(getenv('BENCHMARK_SMOKE_MODE'), '1');
if cfg.smokeMode
    cfg.outputDir = fullfile(figureDir, 'smoke_test_outputs');
else
    cfg.outputDir = fullfile(figureDir, 'outputs');
end
if ~isfolder(cfg.outputDir)
    mkdir(cfg.outputDir);
end

cfg.markovFile = fullfile(cfg.dataDir, 'markov.csv');
cfg.flowFile = fullfile(cfg.dataDir, 'flow.csv');

cfg.nRepeats = 100;
cfg.outerFolds = 5;
cfg.innerFolds = 5;
cfg.bootstrapIterations = 2000;
cfg.randomSeed = 1729;

cfg.lambdaGrid = logspace(-4, 2, 13);
cfg.elasticNetAlpha = 0.5;
cfg.rfMinLeafGrid = [1 3 5 10];
cfg.numTrees = 300;
cfg.primaryTuneMetric = "auc";
cfg.models = ["ridge", "lasso", "elasticnet", "randomforest"];
cfg.views = ["markov", "flow", "clinical_genomic", "markov_plus_clinical"];

if cfg.smokeMode
    cfg.nRepeats = 2;
    cfg.bootstrapIterations = 100;
    cfg.lambdaGrid = [0.1 1 10];
    cfg.rfMinLeafGrid = [1 5];
    cfg.numTrees = 30;
    cfg.models = ["ridge", "lasso"];
    cfg.views = ["markov", "flow"];
end
end
