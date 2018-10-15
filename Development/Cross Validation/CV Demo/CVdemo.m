% CVdemo.m
%
% Training- vs. Test-Set Prediction Error
%
% RTB wrote it, 27 September 2018 procrastinating from T32 study section
% duties

% Basic idea is to generate a data set with a known model, then add noise
% and examine test- vs training-set error as we increase model complexity

%% Generate data using a polynomial of order 3 and add Gaussian noise:

% nPts = 15;
% pTrue = [0, 0.0105, -0.5849, 1.2220];
% x = 1:nPts;
% y = polyval(pTrue,x) + 2*randn(1,nPts);
% yTrue = polyval(pTrue,x);

%% Read in some mystery data

load mysteryData
nPts = length(x);

%% Plot the data:

figure
subplot(3,1,1);
plot(x,y,'bo');
hold on
xlabel('x'); ylabel('y');
title('Mystery Data');

%% Perform fits to data and compare train to test using LOOCV

maxOrder = 7;
trainResid = zeros(nPts,maxOrder);
testResid = zeros(nPts,maxOrder);

warning('off');
for p = 1:maxOrder
    
    % get residuals from fit to all the data:
    pFit = polyfit(x,y,p);
    yPred = polyval(pFit,x);
    trainResid(:,p) = (y - yPred)';
    plot(x,yPred,'-');
    
    % Perform leave-one-out CV to get 'test' residuals
    for k = 1:nPts
        sV = true(1,nPts);      % define a selection vector
        sV(k) = 0;
        
        % fit a polynomial of order p
        pFit = polyfit(x(sV),y(sV),p);
        testResid(k,p) = y(k) - polyval(pFit,x(k));
    end
end
warning('on')

%% Calculate mean squared error from the residuals

MSEtrain = sum(trainResid.^2) ./ nPts;
MSEtest = sum(testResid.^2) ./ nPts;

% figure
subplot(3,1,2);
if (max(MSEtest) - min(MSEtest)) > 10
    semilogy(1:maxOrder,MSEtrain,'b-','LineWidth',2);
    hold on
    semilogy(1:maxOrder,MSEtest,'r-','LineWidth',2);
    legend('Training Sample','Test Sample','Location','Northwest');
else
    plot(1:maxOrder,MSEtrain,'b-','LineWidth',2);
    hold on
    plot(1:maxOrder,MSEtest,'r-','LineWidth',2);
    legend('Training Sample','Test Sample','Location','Northeast');
end

xlabel('Degree of Polynomial');
ylabel('Mean Squared Error');
title('LOOCV')

%% Let's do it with k-fold, where k = 5 or 10

% Repeat CV multiple times to get an estimate of the standard error of the
% estimated test MSE for each model size:
nReps = 50; 
% K-fold cross validation:
kFold = 10;
% Maximum complexity of the polynomial to fit
maxOrder = 7;

% Variables to hold MSEs
MSEtrain = zeros(nReps,maxOrder);   % training error
MSEtest = zeros(nReps,maxOrder);    % cross-validation error

warning('off');
for n = 1:nReps
    for p = 1:maxOrder
        cvIdx = crossvalind('Kfold',nPts,kFold);
        SSEtrain = 0;
        SSEtest = 0;
        for k = 1:kFold
            test = (cvIdx == k);
            train = ~test;
            
            % fit polynomial of order 'p' to training data
            pFit = polyfit(x(train),y(train),p);
            
            % calculate the sum squared error for the training data
            yHatTrain = polyval(pFit,x(train));
            SSEtrain = SSEtrain + sum((yHatTrain - y(train)).^2)/sum(train);
            
            % calculate the sum squared error for the test data
            yHatTest = polyval(pFit,x(test));
            SSEtest = SSEtest + sum((yHatTest - y(test)).^2)/sum(test);
        end
        % Note: This only works if 'kFold' divides evenly into nPts
        % A more general algorithm uses the weighted sum. That is, you
        % calculate the MSE for each fold, then the MSE for each fold
        % contributes to the total in proportion to n_k/n
        MSEtrain(n,p) = SSEtrain / kFold;
        MSEtest(n,p) = SSEtest / kFold;
        
    end
end
warning('on');

% figure
subplot(3,1,3);
if (max(mean(MSEtest)) - min(mean(MSEtest))) > 10
    h1=semilogy(1:maxOrder,mean(MSEtrain),'b-','LineWidth',2);
    hold on
    h2=semilogy(1:maxOrder,mean(MSEtest),'r-','LineWidth',2);
    legFlag = 1;
else
    h1=errorbar(1:maxOrder,mean(MSEtrain),std(MSEtrain),'b-','LineWidth',2);
    hold on
    h2=errorbar(1:maxOrder,mean(MSEtest),std(MSEtest),'r-','LineWidth',2);
    legFlag = 0;
end

xlabel('Degree of Polynomial');
ylabel('Mean Squared Error');
title([num2str(kFold) '-fold Cross Validation']);

%% Choosing the best model

% Note that our value of the test error as a function of the degree of fit
% polynomial suggests that orders of 2, 3 and 4 (and even 5) are rougly
% equivalent in terms of test error. In this setting, we can select a model
% using the one-standard-error rule. We first calculate the standard error
% of the estimated test MSE for each model size, and then select the
% smallest model for which the estimated test error is within one standard
% error of the lowest point on the curve.

% find the point with the lowest test error
minTestError = min(mean(MSEtest));

% draw a horizontal line here:
ax = axis;
hl = line([ax(1),ax(2)],[minTestError, minTestError]);
set(hl,'Color','k','LineStyle','--');

if legFlag
    legend([h1,h2],{'Training Sample','Test Sample'},'Location','Northwest');
else
    legend([h1,h2],{'Training Sample','Test Sample'},'Location','Northeast');
end