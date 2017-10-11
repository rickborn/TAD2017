% etPCAtest.m: Bootstrap calculation of standard error and confidence
% interval for PCA results
%
% RTB adapted for final exam, 11 October 2017

% What to do: Login to learning catalytics and join the session for the
% module entitled "Last Class". You will answer a series of questions based
% on the guided programming below. Each section begins with a '%%'. Read
% through the comments and follow the instructions provided. In some cases
% you will be asked to answer a question, clearly indicated by 'QUESTION'.
% In other cases, you be asked to supply missing code, indicated by 'TODO'.
% The corresponding question in learning catalytics will be indicated in
% parentheses (e.g. Q1). Once you have supplied the required code, you can
% execute that section by mouse-clicking in that section (The block will
% turn yellow.) and then simultaneously hitting the 'ctrl' and 'enter' keys
% (PC) or 'command' and 'enter' keys (Mac). Note that the first section
% doesn't require you to enter any code, so you can just execute this and
% look at the plot. But make sure you read through the comments so that you
% understand the data.
%
% The data here are students' scores on tests in different subjects
% from Mardia, Kent and Bilby (1979)
% n = 88 students (rows) each took 5 tests (columns), in: 
% Mechanics, Vectors, Algebra, Analysis and Statistics. 
% The first two tests were closed book, the last three open book.

% Principal components analysis (PCA) is a useful tool in multivariate
% analysis, as we often have several variables that are correlated with one
% another. PCA allows us to reduce the dimensions of our data by
% transforming our data into orthogonal 'principal components' which
% are uncorrelated and contain the information of the original dataset
% (but may be difficult to interpret). The principal components, 
% or eigenvectors if you've taken linear algebra, are linear combinations 
% of the original variables of the dataset weighted by their contribution 
% to explaining the data's variance in a particular orthogonal dimension,
% chosen such that the first PC explains the most variance possible, the
% second the next most possible in an orthogonal direction to the first,
% and so on.

% Original source of exercise: Efron, B. & Tibshirani Robert, J. (1993) 
% An introduction to the bootstrap. Chapman & Hall, London u.a. pp. 61-70 
% Adapted by RTB 5/31/2002, originally named bs_ex2.m
% Developed for homework by RAS and RTB, August 2017

%% Concepts covered:
%1. z-scoring to see relationships among different measures
%2. PCA to capture maximal variance ('scientific IQ')
%3. bootstrap to estimate std. error of our PCA-based statistic
%4. bootstrap for confidence intervals by the percentile method

%% Load data and make some plots

% NOTE: You will want to make the figure larger so that you can see
% everything clearly.

% Parmeters for bootstrap iterations and significance level. Don't change
% these!
nBoot = 5000;
myAlpha = 0.05;

load MardiaTestScores.mat
X = [Mechanics,Vectors,Algebra,Analysis,Statistics];
[nStudents,nTests] = size(X);

% Plot the raw data to get a sense of the correlations:
figure, subplot(2,2,1);
plot(X);
xlabel('Student ID#'); ylabel('raw test score');
legend('Mechanics','Vectors','Algebra','Analysis','Statistics');
title('Raw test scores');

% We can do even better by first z-scoring the data
% Note that zscore operates on columns by default, so we are essentially
% normalizing each students score within each test
Z = zscore(X);
subplot(2,2,2);
plot(Z);
xlabel('Student ID#'); ylabel('z-scored test score');
legend('Mechanics','Vectors','Algebra','Analysis','Statistics');
title('Z-scored test scores');

% How well correlated are students on two randomly chosen tests
testNames = {'Mechanics','Vectors','Algebra','Analysis','Statistics'};

% Pick 2 random tests, but not the same one twice
testX = unidrnd(nTests);
testY = testX;
while testY == testX
    testY = unidrnd(nTests);
end

subplot(2,2,3);
plot(X(:,testX),X(:,testY),'ko');

xStr = sprintf('Score on test in %s', testNames{testX});
xlabel(xStr);
yStr = sprintf('Score on test in %s', testNames{testY});
ylabel(yStr);
lsline; % plot least squares regression line
title('Test score correlations');

% Look at the entire correlation matrix:
rhoCols = corr(X);
subplot(2,2,4);
imagesc(rhoCols,[0,1]); colorbar;
xlabel('Test'); ylabel('Test');
title('Correlation matrix for test scores');
set(gca, 'XTick', 1:nTests); % center x-axis ticks on bins
set(gca, 'YTick', 1:nTests); % center y-axis ticks on bins
set(gca, 'XTickLabel', testNames);
set(gca, 'YTickLabel', testNames);
set(gca,'XTickLabelRotation',45)

% QUESTION (Q8): Which two tests are most highly correlated with one another?
maxLoc = rhoCols == max(max(tril(rhoCols,-1)));
% Ans. Algebra and analysis

%% PCA to capture sources of maximal variance

% We first note that the test scores are correlated: a student who did well
% on one test also tended to do well on the others. In the limit, each
% student (x) could be captured by a single number, Q, that would
% completely describe her performance: x_i = Q_i * v
%
% We could think of Q as each student's scientific IQ. In the language of
% linear algebra (PCA), v is the eigenvector of the covariance matrix (the
% first principle component) and Q_i is the distance along the line defined
% by the eigenvector for each student. The eigenvalues reflect the amount
% of variance explained by each eigenvector, so the ratio of the 1st
% eigenvalue to the sum of all the eigenvalues gives the proportion of
% variance explained by the first principle component. Thus it tells us how
% much of the variability among the different students can be explained by
% a single number, the scientific IQ.

% Calculate the observed value for the variance explained by the first
% principle component:
% coef: the coefficients for the principal components; 1st PC is column 1, etc.
% score: the representations of X in the principal component space
% latent: the eigenvalues of the covariance matrix: variance explained by each PC
[coef,score,latent] = pca(X);

% TODO: Calculate the proportion of variance explained by the 1st PC and
% store it in 'thHat'
thHat = ;

% QUESTION (Q9): What is the proportion of variance explained by the first
% principle component?

%% Bootstrap to calculate standard error of our measure.

% How precise is our value for proportion of variance explained? There is
% no simple formula for this, like there is for the standard error of the
% mean, but we can do it with the bootstrap by taking random samples of the
% rows, with replacement.

% TODO: Bootstrap the standard error for 'thHat' and store the standard
% error in 'thSE'
% Store your bootstrap replicates in 'thBoot'
thBoot = zeros(nBoot,1);
rng default

!!! Your bootstrap code here

thSE = ;

% QUESTION (Q10): What is your value for the standard error?

% Make a histogram of the bootstrap replicates and draw a vertical line for
% the actual value ('thHat')
figure, hist(thBoot);
hold on
xlabel('Fraction of variance explained by 1st PC');
ylabel('#');
title('Sampling distribution: prop. total var. explained by 1st PC');
ax = axis;
line([thHat,thHat],[ax(3),ax(4)],'Color','y');

%% Bootstrap of confidence intervals: percentile method

% TODO: Use the values in 'thBoot' to calculate a 95% Confidence Interval
% based on the percentile method.
!!! Your code here

% QUESTION (Q11): What is your value for the lower bound of the 95% CI?

% TODO: Draw vertical red lines on your histogram showing the 95% CI
!!! Your code here