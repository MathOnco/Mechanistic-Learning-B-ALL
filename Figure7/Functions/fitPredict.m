function [score, selected, fitted] = fitPredict(modelName, hyperparameter, Xtrain, ytrain, Xtest, cfg)
%FITPREDICT Fit one weighted binary classifier and return P(Y=1).

modelName = lower(string(modelName));
ytrain = double(ytrain(:));
[Xtrain, Xtest] = preprocess(Xtrain, Xtest, modelName ~= "randomforest");
w = classWeights(ytrain);
selected = true(1, size(Xtrain, 2));

switch modelName
    case "ridge"
        Ycat = categorical(ytrain, [0 1], {'Negative','Positive'});
        fitted = fitclinear(Xtrain, Ycat, ...
            'Learner', 'logistic', 'Regularization', 'ridge', ...
            'Lambda', hyperparameter, 'Weights', w, 'Solver', 'lbfgs', ...
            'ClassNames', categorical({'Negative','Positive'}));
        [~, scores] = predict(fitted, Xtest);
        posIdx = fitted.ClassNames == categorical("Positive");
        score = scores(:, posIdx);

    case {"lasso", "elasticnet"}
        if modelName == "lasso"
            alpha = 1;
        else
            alpha = cfg.elasticNetAlpha;
        end
        [B, fitInfo] = lassoglm(Xtrain, ytrain, 'binomial', ...
            'Lambda', hyperparameter, 'Alpha', alpha, 'Weights', w, ...
            'Standardize', false, 'RelTol', 1e-7);
        eta = fitInfo.Intercept + Xtest * B;
        score = 1 ./ (1 + exp(-max(min(eta, 35), -35)));
        if numel(hyperparameter) > 1
            order = zeros(1, numel(hyperparameter));
            used = false(numel(fitInfo.Lambda), 1);
            for i = 1:numel(hyperparameter)
                distance = abs(fitInfo.Lambda(:) - hyperparameter(i));
                distance(used) = Inf;
                [~, k] = min(distance);
                used(k) = true;
                order(i) = k;
            end
            score = score(:, order);
            B = B(:, order);
        end
        selected = abs(B) > 1e-10;
        if isvector(selected)
            selected = selected(:)';
        end
        fitted = struct('Beta', B, 'Intercept', fitInfo.Intercept, 'Alpha', alpha);

    case "randomforest"
        Ycat = categorical(ytrain, [0 1], {'Negative','Positive'});
        fitted = TreeBagger(cfg.numTrees, Xtrain, Ycat, ...
            'Method', 'classification', 'Weights', w, ...
            'MinLeafSize', hyperparameter, 'OOBPrediction', 'off');
        [~, scores] = predict(fitted, Xtest);
        classNames = string(fitted.ClassNames);
        posIdx = classNames == "Positive";
        score = scores(:, posIdx);
        selected = true(1, size(Xtrain, 2));

    otherwise
        error('mechlearn:UnknownModel', 'Unknown classifier: %s', modelName);
end

score = double(score);
if isvector(score)
    score = score(:);
end
score = min(max(score, 0), 1);
end

function [Xtrain, Xtest] = preprocess(Xtrain, Xtest, standardize)
Xtrain = double(Xtrain);
Xtest = double(Xtest);
for j = 1:size(Xtrain, 2)
    finiteTrain = Xtrain(isfinite(Xtrain(:, j)), j);
    if isempty(finiteTrain)
        med = 0;
    else
        med = median(finiteTrain);
    end
    Xtrain(~isfinite(Xtrain(:, j)), j) = med;
    Xtest(~isfinite(Xtest(:, j)), j) = med;
end
if standardize
    mu = mean(Xtrain, 1);
    sigma = std(Xtrain, 0, 1);
    sigma(~isfinite(sigma) | sigma == 0) = 1;
    Xtrain = (Xtrain - mu) ./ sigma;
    Xtest = (Xtest - mu) ./ sigma;
end
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
