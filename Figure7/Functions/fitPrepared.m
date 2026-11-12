function [score, active, modelInfo] = fitPrepared(modelName, lambda, Xtrain, ytrain, Xtest, cfg)
%FITPREPARED Fit one weighted model to already standardized predictors.

modelName = lower(string(modelName));
ytrain = double(ytrain(:));
Xtrain = double(Xtrain);
Xtest = double(Xtest);
w = classWeights(ytrain);

switch modelName
    case "ridge"
        Ycat = categorical(ytrain, [0 1], {'Negative','Positive'});
        fitted = fitclinear(Xtrain, Ycat, ...
            'Learner', 'logistic', 'Regularization', 'ridge', ...
            'Lambda', lambda, 'Weights', w, 'Solver', 'lbfgs', ...
            'ClassNames', categorical({'Negative','Positive'}));
        [~, rawScores] = predict(fitted, Xtest);
        posIdx = string(fitted.ClassNames) == "Positive";
        score = rawScores(:, posIdx);
        beta = double(fitted.Beta(:));
        intercept = double(fitted.Bias);
        active = true(1, size(Xtrain, 2));

    case {"lasso", "elasticnet"}
        if modelName == "lasso", alpha = 1; else, alpha = cfg.elasticNetAlpha; end
        [B, fitInfo] = lassoglm(Xtrain, ytrain, 'binomial', ...
            'Lambda', lambda, 'Alpha', alpha, 'Weights', w, ...
            'Standardize', false, 'RelTol', 1e-7);
        beta = double(B(:, 1));
        intercept = double(fitInfo.Intercept(1));
        eta = intercept + Xtest * beta;
        score = 1 ./ (1 + exp(-max(min(eta, 35), -35)));
        active = abs(beta(:)') > 1e-10;

    otherwise
        error('exhaustive:UnknownModel', 'Unknown model: %s', modelName);
end

score = min(max(double(score(:)), 0), 1);
modelInfo = struct('model', modelName, 'lambda', lambda, ...
    'beta', beta, 'intercept', intercept, 'active', active);
end

function w = classWeights(y)
n = numel(y);
nPos = sum(y == 1);
nNeg = sum(y == 0);
w = ones(n, 1);
if nPos > 0 && nNeg > 0
    w(y == 1) = n / (2 * nPos);
    w(y == 0) = n / (2 * nNeg);
end
end
