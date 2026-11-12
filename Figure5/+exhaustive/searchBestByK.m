function S = searchBestByK(cache, k, p, model, lambda, cfg)
%SEARCHBESTBYK Exhaustively choose the best exact-k candidate in inner CV.

combos = cachedCombinations(p, k);
nCombos = size(combos, 1);
selectionValues = -inf(nCombos, 1);
allMetrics = nan(nCombos, 4);
failed = false(nCombos, 1);
started = tic;

if cfg.useParallel && nCombos > 1
    parfor i = 1:nCombos
        [selectionValues(i), allMetrics(i, :), failed(i)] = ...
            exhaustive.evaluateSubsetInnerCV(cache, combos(i, :), ...
            model, lambda, cfg);
    end
else
    for i = 1:nCombos
        [selectionValues(i), allMetrics(i, :), failed(i)] = ...
            exhaustive.evaluateSubsetInnerCV(cache, combos(i, :), ...
            model, lambda, cfg);
    end
end

bestValue = max(selectionValues);
if ~isfinite(bestValue)
    error('exhaustive:AllSubsetFitsFailed', ...
        'All %d subset fits failed for k=%d.', nCombos, k);
end

ties = find(abs(selectionValues - bestValue) <= cfg.tieTolerance);
if numel(ties) > 1
    tieAuc = allMetrics(ties, 3);
    tieAuc(~isfinite(tieAuc)) = -Inf;
    maxAuc = max(tieAuc);
    ties = ties(abs(tieAuc - maxAuc) <= cfg.tieTolerance);
end
bestIndex = ties(1); % nchoosek order gives deterministic lexicographic tie-break.

S.k = k;
S.cols = combos(bestIndex, :);
S.selectionValue = selectionValues(bestIndex);
S.metricVector = allMetrics(bestIndex, :);
S.nCombinations = nCombos;
S.nFailed = sum(failed);
S.elapsedSeconds = toc(started);
end

function combos = cachedCombinations(p, k)
% Generate each exact-k combination matrix only once per MATLAB session.
persistent cachedP cache
if isempty(cachedP) || cachedP ~= p
    cachedP = p;
    cache = cell(1, p);
end
if isempty(cache{k})
    cache{k} = nchoosek(1:p, k);
end
combos = cache{k};
end
