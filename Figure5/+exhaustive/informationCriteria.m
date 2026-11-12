function I = informationCriteria(modelInfo, Xtrain, ytrain, subsetSize, pTotal, cfg)
%INFORMATIONCRITERIA Descriptive training-set BIC and EBIC values.

ytrain = double(ytrain(:));
eta = modelInfo.intercept + double(Xtrain) * modelInfo.beta;
prob = 1 ./ (1 + exp(-max(min(eta, 35), -35)));
prob = min(max(prob, 1e-12), 1 - 1e-12);
logL = sum(ytrain .* log(prob) + (1 - ytrain) .* log(1 - prob));
n = numel(ytrain);

if modelInfo.model == "ridge"
    Xd = [ones(n, 1), double(Xtrain)];
    workingWeights = max(prob .* (1 - prob), 1e-10);
    H0 = Xd' * (Xd .* workingWeights);
    penalty = eye(size(Xd, 2));
    penalty(1, 1) = 0;
    df = trace((H0 + modelInfo.lambda * penalty) \ H0);
else
    df = 1 + sum(modelInfo.active);
end

bic = -2 * logL + log(n) * df;
logModelCount = gammaln(pTotal + 1) - gammaln(subsetSize + 1) - ...
    gammaln(pTotal - subsetSize + 1);
ebic = bic + 2 .* cfg.ebicGamma .* logModelCount;

I = struct('logLikelihood', logL, 'df', df, 'bic', bic, ...
    'ebicGamma05', ebic(1), 'ebicGamma1', ebic(end));
end
