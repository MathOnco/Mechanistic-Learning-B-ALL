function regenerate_R4_R8_complete_VP(target)
%REGENERATE_R4_R8_COMPLETE_VP Rebuild reviewer Figures R4 and R8.
%
% This correction uses the complete unified virtual-patient table that has:
%   MRD: 104 source observations x 100 fits = 10,400 rows
%   BCR: 101 source observations x 100 fits = 10,100 rows
%
% Figure R8 assigns folds by source specimen: every selected realization
% of a specimen stays entirely in training or entirely in test. Results
% are written to outputs/ inside this folder.
%
% Usage:
%   regenerate_R4_R8_complete_VP          % build both figures
%   regenerate_R4_R8_complete_VP("R4")    % build R4 only
%   regenerate_R4_R8_complete_VP("R8")    % build/resume R8 only

if nargin < 1
    target = "both";
end
target = upper(string(target));
if ~ismember(target, ["BOTH", "R4", "R8"])
    error('Target must be "both", "R4", or "R8".');
end

thisDir = fileparts(mfilename('fullpath'));
completeVPFile = fullfile(thisDir, 'Data', 'PB_BM_merged_with_MRDflow_BCR.csv');
outputDir = fullfile(thisDir, 'outputs');
figureDir = fullfile(outputDir, 'response_figures');

assert(isfile(completeVPFile), 'Complete VP table not found: %s', completeVPFile);
if ~isfolder(outputDir), mkdir(outputDir); end
if ~isfolder(figureDir), mkdir(figureDir); end
addpath(thisDir);
addpath(fullfile(thisDir, 'Functions'));

V = loadAndValidateCompleteVP(completeVPFile);
fprintf(['Validated complete VP table: %d blocks, %d rows; ' ...
    'MRD=%d x 100, BCR=%d x 100.\n'], ...
    V.nBlocks, height(V.T), V.nMRDBlocks, V.nBCRBlocks);

P.sourceDir = outputDir;
P.figureDir = figureDir;
P.isSmoke = false;
P.sourceLabel = "complete virtual-patient analysis";

if target == "BOTH" || target == "R4"
    makeCorrectedR4(V, outputDir, P);
end
if target == "BOTH" || target == "R8"
    makeCorrectedR8(V, outputDir, P);
end

fprintf('\nCorrected outputs are in:\n%s\n', outputDir);
end

function V = loadAndValidateCompleteVP(filename)
T = readtable(filename, 'VariableNamingRule', 'preserve', 'TextType', 'string');

featureNames = ["M11","M12","M13","M14","M21","M22","M23","M24", ...
    "M31","M32","M33","M34","M41","M42","M43","M44"];
required = ["ID","VirtualPatient",featureNames,"x1","x3", ...
    "MRD(flow)","BCR::ABL1/like"];
missing = setdiff(required, string(T.Properties.VariableNames));
if ~isempty(missing)
    error('Complete VP table is missing: %s', strjoin(missing, ', '));
end

id = double(T.ID);
vp = double(T.VirtualPatient);
if isempty(vp) || vp(1) ~= 1
    error('The first virtual-patient row must have VirtualPatient=1.');
end
newBlock = vp == 1;
block = cumsum(newBlock);
nBlocks = max(block);

blockSizes = accumarray(block, 1, [nBlocks 1]);
if any(blockSizes ~= 100)
    bad = find(blockSizes ~= 100);
    error('Every source block must contain 100 rows. Invalid block(s): %s', ...
        strjoin(string(bad'), ', '));
end
for b = 1:nBlocks
    rows = block == b;
    expected = (1:100)';
    if ~isequal(vp(rows), expected)
        error('Block %d does not contain VirtualPatient values 1:100 in order.', b);
    end
    if numel(unique(id(rows))) ~= 1
        error('Block %d contains more than one source ID.', b);
    end
end

mrd = strtrim(string(T{:, 'MRD(flow)'}));
bcr = strtrim(string(T{:, 'BCR::ABL1/like'}));
mrdEligible = ismember(mrd, ["Negative","Positive (< 10^-3)", ...
    "Positive (> 10^-3)"]);
bcrEligible = ismember(bcr, ["Negative","BCR::ABL1"]);

mrdBlocks = unique(block(mrdEligible), 'stable');
bcrBlocks = unique(block(bcrEligible), 'stable');
if numel(mrdBlocks) ~= 104 || sum(mrdEligible) ~= 10400
    error('Expected 104 MRD blocks/10,400 rows; found %d/%d.', ...
        numel(mrdBlocks), sum(mrdEligible));
end
if sum(mrd == "Positive (> 10^-3)") ~= 1900 || ...
        sum(mrd == "Negative" | mrd == "Positive (< 10^-3)") ~= 8500
    error('MRD class counts do not equal 1,900 positive and 8,500 negative.');
end
if numel(bcrBlocks) ~= 101 || sum(bcrEligible) ~= 10100
    error('Expected 101 BCR blocks/10,100 rows; found %d/%d.', ...
        numel(bcrBlocks), sum(bcrEligible));
end
if sum(bcr == "BCR::ABL1") ~= 5100 || sum(bcr == "Negative") ~= 5000
    error('BCR class counts do not equal 5,100 positive and 5,000 negative.');
end

% Confirm that each outcome label is constant within a 100-fit source block.
for b = 1:nBlocks
    rows = block == b;
    m = unique(mrd(rows));
    m = m(~ismissing(m) & strlength(strtrim(m)) > 0);
    c = unique(bcr(rows));
    c = c(~ismissing(c) & strlength(strtrim(c)) > 0);
    if numel(m) > 1 || numel(c) > 1
        error('Outcome label changes within source block %d.', b);
    end
end

V.T = T;
V.id = id;
V.vp = vp;
V.block = block;
V.nBlocks = nBlocks;
V.blockSizes = blockSizes;
V.featureNames = featureNames;
V.X = double(T{:, cellstr(featureNames)});
V.x1 = double(T.x1);
V.x3 = double(T.x3);
V.mrd = mrd;
V.bcr = bcr;
V.mrdEligible = mrdEligible;
V.bcrEligible = bcrEligible;
V.nMRDBlocks = numel(mrdBlocks);
V.nBCRBlocks = numel(bcrBlocks);
end

function makeCorrectedR4(V, outputDir, P)
n = V.nBlocks;
sourceSample = (1:n)';
patientID = nan(n, 1);
nFits = nan(n, 1);
mrdFlow = strings(n, 1);
bcrStatus = strings(n, 1);
meanM11 = nan(n, 1); sdM11 = nan(n, 1);
meanM14 = nan(n, 1); sdM14 = nan(n, 1);
meanM31 = nan(n, 1); sdM31 = nan(n, 1);
sourceAbundanceX1 = nan(n, 1);
sourceAbundanceX3 = nan(n, 1);
rangeX1 = nan(n, 1); rangeX3 = nan(n, 1);

iM11 = find(V.featureNames == "M11", 1);
iM14 = find(V.featureNames == "M14", 1);
iM31 = find(V.featureNames == "M31", 1);
for b = 1:n
    rows = V.block == b;
    first = find(rows, 1);
    patientID(b) = V.id(first);
    nFits(b) = sum(rows);
    mrdFlow(b) = firstNonempty(V.mrd(rows));
    bcrStatus(b) = firstNonempty(V.bcr(rows));
    meanM11(b) = mean(V.X(rows, iM11), 'omitnan');
    sdM11(b) = std(V.X(rows, iM11), 0, 'omitnan');
    meanM14(b) = mean(V.X(rows, iM14), 'omitnan');
    sdM14(b) = std(V.X(rows, iM14), 0, 'omitnan');
    meanM31(b) = mean(V.X(rows, iM31), 'omitnan');
    sdM31(b) = std(V.X(rows, iM31), 0, 'omitnan');
    sourceAbundanceX1(b) = mean(V.x1(rows), 'omitnan');
    sourceAbundanceX3(b) = mean(V.x3(rows), 'omitnan');
    rangeX1(b) = max(V.x1(rows)) - min(V.x1(rows));
    rangeX3(b) = max(V.x3(rows)) - min(V.x3(rows));
end
if max(rangeX1, [], 'omitnan') > 1e-12 || max(rangeX3, [], 'omitnan') > 1e-12
    warning('Source abundances vary within at least one VP block; block means were used.');
end

S = table(sourceSample, patientID, nFits, mrdFlow, bcrStatus, ...
    meanM11, sdM11, meanM14, sdM14, meanM31, sdM31, ...
    sourceAbundanceX1, sourceAbundanceX3);
writetable(S, fullfile(outputDir, 'transition_variability_by_source_sample.csv'));

transition = ["M11";"M14";"M31"];
nValid = zeros(3, 1); spearmanRho = nan(3, 1); pValue = nan(3, 1);
abundance = {sourceAbundanceX1, sourceAbundanceX1, sourceAbundanceX3};
uncertainty = {sdM11, sdM14, sdM31};
for j = 1:3
    valid = isfinite(abundance{j}) & isfinite(uncertainty{j});
    nValid(j) = sum(valid);
    [spearmanRho(j), pValue(j)] = corr(abundance{j}(valid), ...
        uncertainty{j}(valid), 'Type', 'Spearman');
end
C = table(transition, nValid, spearmanRho, pValue, ...
    'VariableNames', {'transition','n','spearmanRho','pValue'});
writetable(C, fullfile(outputDir, 'source_abundance_vs_transition_SD.csv'));

figureTransitionUncertainty(P);
fprintf('R4 regenerated from %d complete 100-fit source blocks.\n', n);
end

function makeCorrectedR8(V, outputDir, P)
cfg = reviewer_config();
% Force publication settings even if REVIEWER_SMOKE_MODE is set elsewhere.
cfg.smokeMode = false;
cfg.outerFolds = 5;
cfg.innerFolds = 5;
cfg.randomSeed = 1729;
cfg.lambdaGrid = logspace(-4, 2, 13);
cfg.elasticNetAlpha = 0.5;
cfg.primaryTuneMetric = "auc";
cfg.models = ["ridge","lasso"];
cfg.vpCounts = [1 5 10 25 50 100];
cfg.vpRegimeRepeats = 30;
cfg.vectorizedLambda = true;
warning('off', 'stats:lasso:MaxIterReached');
warning('off', 'stats:glmfit:IterationLimit');

metricsFile = fullfile(outputDir, ...
    'VP_augmentation_regime_metrics_specimen_folds.csv');
if isfile(metricsFile)
    allMetrics = readtable(metricsFile, 'TextType', 'string', ...
        'VariableNamingRule', 'preserve');
    fprintf('Resuming R8 from %s (%d saved rows).\n', metricsFile, height(allMetrics));
else
    allMetrics = table();
end

functionsDir = fileparts(which('nestedCV'));
if license('test', 'Distrib_Computing_Toolbox')
    pool = gcp('nocreate');
    if isempty(pool)
        pool = parpool;
    end
    pctRunOnAll(['addpath(''' strrep(functionsDir, '''', '''''') ''')']);
    fprintf('Parallel pool: %d workers\n', pool.NumWorkers);
end

outcomes = ["MRDflow","BCR"];
for o = 1:numel(outcomes)
    outcome = outcomes(o);
    if outcome == "MRDflow"
        eligible = V.mrdEligible;
        yAll = double(V.mrd == "Positive (> 10^-3)");
        expectedBlocks = 104;
    else
        eligible = V.bcrEligible;
        yAll = double(V.bcr == "BCR::ABL1");
        expectedBlocks = 101;
    end
    eligibleBlocks = unique(V.block(eligible), 'stable');
    assert(numel(eligibleBlocks) == expectedBlocks);

    for q = cfg.vpCounts
        pending = 1:cfg.vpRegimeRepeats;
        if ~isempty(allMetrics)
            for r = pending
                done = allMetrics.outcome == outcome & ...
                    allMetrics.virtualPatientsPerSource == q & ...
                    allMetrics.samplingReplicate == r;
                if sum(done) == numel(cfg.models)
                    pending(pending == r) = [];
                elseif any(done)
                    allMetrics(done, :) = [];
                end
            end
        end
        if isempty(pending)
            continue;
        end
        batch = cell(numel(pending), 1);
        parfor ii = 1:numel(pending)
            batch{ii} = runVPReplicate(V, yAll, eligible, eligibleBlocks, ...
                outcome, q, pending(ii), o, cfg);
        end
        for ii = 1:numel(batch)
            allMetrics = [allMetrics; batch{ii}]; %#ok<AGROW>
        end
        allMetrics = sortrows(allMetrics, ...
            {'outcome','virtualPatientsPerSource','samplingReplicate','model'});
        writetable(allMetrics, metricsFile);
        fprintf('%s q=%3d finished %d new replicates.\n', ...
            outcome, q, numel(pending));
    end
end

[G, outcomeG, qG, modelG] = findgroups(allMetrics.outcome, ...
    allMetrics.virtualPatientsPerSource, allMetrics.model);
summary = table(outcomeG, qG, modelG, ...
    splitapply(@mean, allMetrics.auc, G), ...
    splitapply(@(x) prctile(x, 2.5), allMetrics.auc, G), ...
    splitapply(@(x) prctile(x, 97.5), allMetrics.auc, G), ...
    splitapply(@mean, allMetrics.balancedAccuracy, G), ...
    splitapply(@(x) prctile(x, 2.5), allMetrics.balancedAccuracy, G), ...
    splitapply(@(x) prctile(x, 97.5), allMetrics.balancedAccuracy, G), ...
    'VariableNames', {'outcome','virtualPatientsPerSource','model', ...
    'meanAUC','aucLower','aucUpper','meanBalancedAccuracy', ...
    'balaccLower','balaccUpper'});
writetable(summary, fullfile(outputDir, ...
    'VP_augmentation_regime_summary_specimen_folds.csv'));

figureVPAugmentation(P);
fprintf('R8 regenerated with 104 MRD and 101 BCR source blocks.\n');
end

function M = runVPReplicate(V, yAll, eligible, eligibleBlocks, ...
        outcome, q, replicate, outcomeIndex, cfg)
warning('off', 'stats:lasso:MaxIterReached');
warning('off', 'stats:glmfit:IterationLimit');
rng(cfg.randomSeed + 100000 * outcomeIndex + 1000 * q + replicate, 'twister');
selectedRows = zeros(numel(eligibleBlocks) * q, 1);
cursor = 0;
for b = eligibleBlocks'
    candidates = find(V.block == b & eligible);
    if numel(candidates) ~= 100
        error('Eligible source block %d does not contain 100 rows.', b);
    end
    chosen = candidates(randperm(100, q));
    selectedRows(cursor + (1:q)) = chosen;
    cursor = cursor + q;
end
selectedRows = selectedRows(randperm(numel(selectedRows)));

localCfg = cfg;
localCfg.nRepeats = 1;
label = outcome + "__VP_count_" + q + "__sample_" + replicate;
R = nestedCV(V.X(selectedRows, :), yAll(selectedRows), ...
    V.featureNames, cfg.models, localCfg, label, V.block(selectedRows));
M = R.repeatMetrics;
M.outcome = repmat(outcome, height(M), 1);
M.virtualPatientsPerSource = repmat(q, height(M), 1);
M.samplingReplicate = repmat(replicate, height(M), 1);
M.nSourceSamples = repmat(numel(eligibleBlocks), height(M), 1);
M.splitUnit = repmat("source_specimen", height(M), 1);
end

function value = firstNonempty(values)
values = string(values(:));
values = values(~ismissing(values) & strlength(strtrim(values)) > 0);
if isempty(values)
    value = "";
else
    value = values(1);
end
end
