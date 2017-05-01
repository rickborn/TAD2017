% etPCAdemo.m: Bootstrap calculation of std. error for PCA results
%
% from pp. 61-70 of Efron & Tibshirani
% RTB 5/31/2002, originally named bs_ex2.m
%
% Demo version of etPCAstats.m (function version)

% Concepts covered:
1. z-scoring to see relationships among different measures
2. PCA to capture maximal variance ('scientific IQ')
3. bootstrap to estimate std. error of our PCA-based statistic
4. bootstrap for confidence intervals: std. normal vs. percentile

%% Load and plot data
nBoot = 200;
myAlpha = 0.05;

% Test score data from Mardia, Kent and Bilby (1979)
% n = 88 students (rows) each took 5 tests (columns), in: 
% Mechanics, Vectors, Algebra, Analysis and Statistics. 
% The first two tests were closed book, the last three open book.
load MardiaTestScores.mat
X = [Mechanics,Vectors,Algebra,Analysis,Statistics];
[nStudents,nTests] = size(X);

% Plot the raw data to get a sense of the correlations:
figure, plot(X);
xlabel('Student ID#'); ylabel('raw test score');
legend('Mechanics','Vectors','Algebra','Analysis','Statistics');
title('Raw test scores from Mardia et al. 1979');

% We can do even better by first z-scoring the data
% Note that zscore operates on columns by default, so we are essentially
% normalizing each students score within each test
Z = zscore(X);
figure, plot(Z);
xlabel('Student ID#'); ylabel('z-scored test score');
legend('Mechanics','Vectors','Algebra','Analysis','Statistics');
title('Z-scored test scores from Mardia et al. 1979');
%% How well correlated are students on two randomly chosen tests
testNames = {'Mechanics','Vectors','Algebra','Analysis','Statistics'};
testX = unidrnd(nTests);
testY = unidrnd(nTests);
figure, plot(X(:,testX),X(:,testY),'ko');
xStr = sprintf('Score on test in %s', testNames{testX});
xlabel(xStr);
yStr = sprintf('Score on test in %s', testNames{testY});
ylabel(yStr);
lsline; % plot least squares regression line
title('Correlation between two randomly selected tests');
%% Also look at correlation matrix:
rhoCols = corr(X);
figure, imagesc(rhoCols,[0,1]); colorbar;
xlabel('Test #'); ylabel('Test #');
title('Correlation matrix for test scores');

%% PCA to capture sources of maximal variance
% We first note that the test scores are highly correlated: a student who
% did well on one test also tended to well on the others. In the limit,
% each student could be captured by a single number, Q, that would
% completely describe his performance: x_i = Q_i * v
% We could think of Q as each student's scientific IQ.
% In the language of linear algebra (PCA), v is the eigenvector of the
% covariance matrix (the first principle component) and Q_i is the distance
% along the line defined by the eigenvector for each student. The
% eigenvalues reflect the amount of variance explained by each eigenvector,
% so the ratio of the 1st eigenvalue to the sum of all the eigenvalues
% gives the fraction of variance explained by the first principle
% component. Thus it tells us how much of the variability among the
% different stuudents can be explained by a single number, the IQ.

% Calculate the observed value for the variance explained by the first
% principle component:
% coef: the coeficients for the principal components; 1st PC is column 1, etc.
% score: the representations of X in the principal component space
% latent: the eigenvalues of the covariance matrix; related to variance
% explained by each PD
[coef,score,latent] = pca(X);
% thHat is the proportion of the total variance explained by the 1st PC
thHat = latent(1) / sum(latent);    % 0.6191

%% Bootstrap to calculate standard error of our measure.
% So the first PC explains 62% of the test score variance, which is pretty
% darned good. But how precise is this value? There is no simple formula
% for this, like there is for the std. error of the mean, but we can do it
% with the bootstrap by taking random samples of the rows, with
% replacement:
thBoot = zeros(nBoot,1);
coefBoot = zeros(nTests,nTests,nBoot);
for k = 1:nBoot
    bootSamp = X(unidrnd(nStudents,nStudents,1),:);
    [bsCoef,~,bsLatent] = pca(bootSamp);
    thBoot(k) = bsLatent(1) / sum(bsLatent);
    % also save the coeficients for the PCs
    coefBoot(:,:,k) = bsCoef;
end
figure, hist(thBoot);
hold on
xlabel('Fraction of variance explained by 1st PC');
ylabel('#');
title('Bootstrap replicates of proportion of the total variance explained by the 1st PC');
thSE = std(thBoot);
ax = axis;
line([thHat,thHat],[ax(3),ax(4)],'Color','y');

%% Bootstrap of confidence intervals
% There are two ways to compute confidence intervals with the bootstrap.
% One is to do a huge number of bootstrap replicates and then sort them and
% find the cut-offs. E.g. If we do 1000 replicates, and we want the 95% CI,
% we sort and find the 25th and 975th values.
% The other way is to use asymptotic normal distribution theory. Recall
% that the std. error IS a confidence interval, the 68% CI. So to get any
% other CI, we just calculate the appropriate number of standard deviates
% from the normal distribution:
dev = thSE * norminv(1-myAlpha/2);
thCI = [thHat - dev, thHat + dev];
h1 = line([thCI(1),thCI(1)],[ax(3),ax(4)],'Color','r');
line([thCI(2),thCI(2)],[ax(3),ax(4)],'Color','r');

%% Here is the brute force way: percentile method
thBootSorted = sort(thBoot);
iLo = ceil((myAlpha/2) * nBoot);   % index corresponding to lower bound
iHi = nBoot - iLo;                  % index corresponding to upper bound
thCI2 = [thBootSorted(iLo),thBootSorted(iHi)];
h2=line([thCI2(1),thCI2(1)],[ax(3),ax(4)],'Color','g');
line([thCI2(2),thCI2(2)],[ax(3),ax(4)],'Color','g');

legend([h1,h2],'Standard Normal CI','Percentile CI');

%% We can also look at the variability of the PC coeficients themselves:
figure,
h1=plot(squeeze(coefBoot(:,1,:)),'k-');         % 1st PC
hold on
h2=plot(squeeze(coefBoot(:,2,:)),'r-');            % 2nd PC
h3=plot(squeeze(coefBoot(:,3,:)),'g-');            % 3rd PC
legend([h1(1),h2(1),h3(1)],'PC1','PC2','PC3');
xlabel('Coefficient #'); ylabel('Coefficient value');
title('Bootstrap replicates of 1st three PCs');

%% Computation of the PC's.
% Start with the covariance matrix:
G = cov(X);
% The PC are obtained as the eigenvectors of the covariance matrix:
[V,D] = eig(G);

% The eigenvalues give a relative measure of the variance explained by each
% PC. In MATLAB's PCA-speak, this is returned as 'latent'. But we can also
% get it as:
lambda = diag(D);   % lambda for the 1st PC is lambda(5), seemingly backwards

% For E&T p. 67-8 on interpretation of PCs:
% First eigenvector (PC1) = V(:,end)' = [0.505,0.368,0.3456,0.451,0.535]
% Second eigenvector (PC2) = V(:,end-1)' = [-0.75,-0.21,0.08,0.30,0.55]
%
% If we wanted to best summarize each student's performance with one
% number, we would multiply each test score by the corresponding
% coefficient and sum up the result. This is called the 'loading':
loadPC1 = X(:,:) * V(:,5);  % Scientific IQ for each of the 88 students
loadPC2 = X(:,:) * V(:,4);  % loading onto 2nd PC

% The weights assigned by PCA can tell us about the structure of the
% multivariate data set. The 1st PC puts positive and roughly equal weights
% on each of the five tests, so 'loadPC1' can be thought of as the sum (or
% average) of each student's 5 test scores. PC2 puts negative weights on
% the first 2 (closed book) tests and positive weights on the 3 open-book
% tests, so loadPC2 is a contrast between a student's open and close book
% performances. A student with a high 'loadPC2' did much better on the open
% book tests than on the closed book tests.

