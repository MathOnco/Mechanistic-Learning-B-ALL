function [Xtrain, Xtest, prep] = prepareData(Xtrain, Xtest)
%PREPAREDATA Training-only median imputation and standardization.

Xtrain = double(Xtrain);
Xtest = double(Xtest);
p = size(Xtrain, 2);
med = zeros(1, p);
for j = 1:p
    finiteTrain = Xtrain(isfinite(Xtrain(:, j)), j);
    if isempty(finiteTrain), med(j) = 0; else, med(j) = median(finiteTrain); end
    Xtrain(~isfinite(Xtrain(:, j)), j) = med(j);
    Xtest(~isfinite(Xtest(:, j)), j) = med(j);
end
mu = mean(Xtrain, 1);
sigma = std(Xtrain, 0, 1);
sigma(~isfinite(sigma) | sigma == 0) = 1;
Xtrain = (Xtrain - mu) ./ sigma;
Xtest = (Xtest - mu) ./ sigma;
prep = struct('median', med, 'mean', mu, 'sd', sigma);
end
