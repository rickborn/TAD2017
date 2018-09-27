% CVfail.m: example of the wrong way to use cross-validation
%
% RTB wrote it, 06 July 2018

% Based on an example from the Hastie-Tibshirani MOOC on statistical
% learning:
%
% Imagine we want a simple classifier for some 2-class data: 
% 1) Starting with 5000 predictors and 50 samples, we first filter by 
%    finding the 100 predictors having the largest correlation with the 
%    class labels 
% 2) We then apply a classifier, such as logistic regression, using only 
%    these 100 predictors.
% 
% Can we apply CV to step #2, forgetting step #1?

%% Simulate some data

nPredictors = 5000;     % # of predictor columns in raw data
nSamples = 100;         % # of samples
nCherry = 50;           % # of best predictors to cherry pick

% Note: To use LDA without regularization, the number of predictors must be
% less than the number of samples (taking into acount the reduction in the
% training set size created by cross-validation). If you violate this rule,
% 'classify' will return the following error message:
%
% Error using classify (line 233)
% The pooled covariance matrix of TRAINING must be positive definite.

% Generate random data, in which there is no true predictive value
D = randn(nSamples,nPredictors);
% Arbitrarily assign half of the samples to class 0; half to class 1
dLabels = [zeros(nSamples/2,1);ones(nSamples/2,1)];

%% CV the wrong way

% Step 1: filter by finding the best predictors of class
allR = zeros(nPredictors,1);

% The clever matrix way to do this is slower than the 'for' loop, because
% you are computing correlations among all columns. But here is how you
% would do it:
% E = [dLabels,D];
% allRfull = corrcoef(E);
% allR2 = allRfull(2:end,1);
for k = 1:nPredictors
    allR(k,1) = diag(corrcoef(dLabels,D(:,k)),1);
end

% Now find the nCherry best predictors:
[~,I] = sort(abs(allR),'descend');
bestD = D(:,I(1:nCherry,1));

%% Do K-fold cross-validation
K = 10;
indices = crossvalind('Kfold', nSamples, K);
cp = classperf(dLabels);

% Note: We could use a variety of methods to classify the data. For
% example, linear or logistic regression would work. In this example, I'm
% using MATLAB's 'classify' function. It's default mode is to do linear
% discriminant analysis. But it can do other varieties as well. For
% quadratic discriminant use 'quadratic'; for naive Bayes, use
% 'diagquadratic'
for k = 1:K
    test = (indices == k);
    train = ~test;
    class = classify(bestD(test,:),bestD(train,:),dLabels(train,:),'linear');
    classperf(cp,class,test);
end

% cp.ErrorRate is the fraction of test samples classified incorrectly
fprintf(1,'The fraction of test samples incorrectly classified is %0.2f.\n',...
    cp.ErrorRate);

% This gives us a very low error rat (< 15%) even though there is no useful
% information in our predictors. Why? Because we have cheated with the
% pre-filtering and selection of predictors based on ALL of the data (i.e.
% before the CV step.

% What should the error rate be?

%% CV the right way

% The correct way to do this is to partition the data set prior to
% picking the best predictors. That way, the best predictors are chosen
% based only on the training data, and this will provide no useful
% prediction for the test data.

% Do K-fold cross-validation
K = 10;
indices = crossvalind('Kfold', nSamples, K);
cp = classperf(dLabels);

for k = 1:K
    test = (indices == k);
    train = ~test;
    
    % filter, but only use the training data
    allR = zeros(nPredictors,1);
    for k = 1:nPredictors
        allR(k,1) = diag(corrcoef(dLabels(train,1),D(train,k)),1);
    end
    
    % Now find the best predictors:
    [~,I] = sort(abs(allR),'descend');
    bestD = D(:,I(1:nCherry,1));
    
    % Finally, run the classifier and test performance
    class = classify(bestD(test,:),bestD(train,:),dLabels(train,:));
    classperf(cp,class,test);
end

% cp.ErrorRate is the fraction of test samples classified incorrectly
fprintf(1,'The fraction of test samples incorrectly classified is %0.2f.\n',...
    cp.ErrorRate);

% Now we get error rates around 50%, which is what we should get!