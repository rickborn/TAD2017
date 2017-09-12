% etPCAdemo.m: Bootstrap calculation of std. error for PCA results
%
% Instructions:
% The goal of this exercise is to become familiar with the technique of
% bootstrapping and appreciate how it can be used to estimate accuracy of 
% statistics through resampling data to generate standard errors and
% confidence intervals that may otherwise be difficult to compute directly.

% What to do: Login to learning catalytics and join the session for the
% module entitled "etPCAdemo, self-test". You will answer a series of
% questions based on the guided programming below. Each section begins with
% a '%%'. Read through the comments and follow the instructions provided.
% In some cases you will be asked to answer a question, clearly indicated
% by 'QUESTION'. In other cases, you be asked to supply missing code,
% indicated by 'TODO'. The corresponding question in learning catalytics
% will be indicated in parentheses (e.g. Q1). If there is no 'Q#'
% accompanying a 'QUESTION' just type your answer into this script and
% discuss it with your team. Once you have supplied the required code, you
% can execute that section by mouse-clicking in that section (The block
% will turn yellow.) and then simultaneously hitting the 'ctrl' and 'enter'
% keys (PC) or 'command' and 'enter' keys (Mac).
%
% The data here are students' scores across tests in different subjects
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
%4. bootstrap for confidence intervals: std. normal vs. percentile

%% Load and plot data
nBoot = 5000;
myAlpha = 0.05;

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

% QUESTION (Q1): Looking at Figure 1 and Figure 2, does there seem to be any
% relationship between a student's scores on the different tests? How would
% you get a sense for how well different tests track together on this type
% of plot?
%% How well correlated are students on two randomly chosen tests

testNames = {'Mechanics','Vectors','Algebra','Analysis','Statistics'};
% TO DO: We want to randomly choose two tests to plot against each other.
% Define testX and testY to be variables that will be a random whole number
% from 1 to the number of tests (this information is stored in the variable
% nTests).
testX = 
testY = 

% TO DO: Write a single line of code to plot the correlation between two
% randomly selected tests. Use testX and testY to index into the main data
% array X. Make the points black and open circles. The lines after will set
% up axes and titles.
figure
% YOUR CODE HERE

xStr = sprintf('Score on test in %s', testNames{testX});
xlabel(xStr);
yStr = sprintf('Score on test in %s', testNames{testY});
ylabel(yStr);
lsline; % plot least squares regression line
title('Correlation between two randomly selected tests');

%QUESTION (Q2): Run this block of code several times. Do tests appear to be
%correlated (again this is only by eye at this point). In what direction?
%% Also look at correlation matrix:

% TODO: Use the 'corr' function to create a matrix that shows the
% correlation coefficient for each test with every other. Use 'imagesc' to
% diplay the matrix, using the color to indicate the correlation
% coefficient between each pair of tests:
rhoCols = ;
figure
% YOUR CODE HERE

% Label the plot
xlabel('Test'); ylabel('Test');
title('Correlation matrix for test scores');
set(gca, 'XTick', 1:nTests); % center x-axis ticks on bins
set(gca, 'YTick', 1:nTests); % center y-axis ticks on bins
set(gca, 'XTickLabel', testNames);
set(gca, 'YTickLabel', testNames);

% QUESTION (Q3): Which tests appear to be the most correlated with one another?
% In general, are the test scores highly correlated with one another?

%% PCA to capture sources of maximal variance
% We first note that the test scores are correlated: a student who
% did well on one test also tended to do well on the others. In the limit,
% each student (x) could be captured by a single number, Q, that would
% completely describe their performance: x_i = Q_i * v
%
% We could think of Q as each student's scientific IQ.
% In the language of linear algebra (PCA), v is the eigenvector of the
% covariance matrix (the first principle component) and Q_i is the distance
% along the line defined by the eigenvector for each student. The
% eigenvalues reflect the amount of variance explained by each eigenvector,
% so the ratio of the 1st eigenvalue to the sum of all the eigenvalues
% gives the fraction of variance explained by the first principle
% component. Thus it tells us how much of the variability among the
% different students can be explained by a single number, the IQ.

% Calculate the observed value for the variance explained by the first
% principle component:
% coef: the coefficients for the principal components; 1st PC is column 1, etc.
% score: the representations of X in the principal component space
% latent: the eigenvalues of the covariance matrix; related to variance
% explained by each PC
[coef,score,latent] = pca(X);

% TO DO: Using the knowledge above, assign the variable thHat 
% to be the proportion of the total variance explained by the 1st PC.
thHat = ;

% QUESTION (Q4): What is the proportion of the total variance explained by
% the 1st PC?

%% Bootstrap to calculate standard error of our measure.
% How precise is the value we calculated above? There is no simple formula
% for this, like there is for the standard error of the mean, but we can do
% it with the bootstrap by taking random samples of the rows, with
% replacement:

% Initialize arrays for holding the bootstrap replicates:
thBoot = zeros(nBoot,1);
coefBoot = zeros(nTests,nTests,nBoot);

rng default
for k = 1:nBoot
    % TO DO: write a single line of code to randomly sample the rows of our
    % dataset X with replacement. The resulting variable bootSamp, will be
    % an 88x5 double.
    % HINT: This requires careful indexing of our dataset X and can be 
    % accomplished with the use of the unidrnd function to draw numbers
    % from a random uniform distribution with the maximum value set at the
    % number of students.
    bootSamp = ;
    
    [bsCoef,~,bsLatent] = pca(bootSamp);
    thBoot(k) = bsLatent(1) / sum(bsLatent);
    % also save the coefficients for the PCs
    coefBoot(:,:,k) = bsCoef;
end
figure, hist(thBoot);
hold on
xlabel('Fraction of variance explained by 1st PC');
ylabel('#');
title('Bootstrap replicates of proportion of the total variance explained by the 1st PC');

% Recall that the standard error of any statistic is defined as the
% standard deviation of its sampling distribution:
thSE = std(thBoot);

% Line showing the value we actually obtained:
ax = axis;
line([thHat,thHat],[ax(3),ax(4)],'Color','y');

%% Bootstrap of confidence intervals
% There are two general approaches to computing confidence intervals with the bootstrap:
%
% One way is to use asymptotic normal distribution theory. Recall
% that the std. error IS a confidence interval, the 68% CI. So to get any
% other CI, we just calculate the appropriate number of standard deviates
% from the normal distribution:
dev = thSE * norminv(1-myAlpha/2);
thCI = [thHat - dev, thHat + dev];
h1 = line([thCI(1),thCI(1)],[ax(3),ax(4)],'Color','r');
line([thCI(2),thCI(2)],[ax(3),ax(4)],'Color','r');

% QUESTION (Q5): How would you modify the above code if you wanted to get a
% 99% CI?

%% CI by percentile method
% The other way is to do a huge number of bootstrap replicates and then
% sort them and find the cut-offs. E.g. If we do 1000 replicates, and we
% want the 95% CI, we sort and find the 25th and 975th values.

% TODO: Using your bootstrap replicates in 'thBoot', use the percential
% method to calculate the lower and upper bounds for the 95% CI. Store them
% in 'thCIlow' and 'thCIhi', respectively:
thCIlow = ;
thCIhi = ;
thCI2 = [thBootSorted(iLo),thBootSorted(iHi)];
h2=line([thCIlow, thCIlow],[ax(3),ax(4)],'Color','g');
line([thCIhi, thChi],[ax(3),ax(4)],'Color','g');

% QUESTION (Q6): What is the value of your lower bound ('thCIlow')?

legend([h1,h2],'Standard Normal CI','Percentile CI');
% QUESTION (Q7): How similar are the standard normal CI and percentile CIs? 
% Why might they be slightly different? 

%% We can also look at the variability of the PC coefficients themselves:

% subsample of 200 for viewing clarity:
whichRows = randi(nBoot,200,1);

figure,
h1=plot(squeeze(coefBoot(:,1,whichRows)),'k-'); % 1st PC
hold on
h2=plot(squeeze(coefBoot(:,2,whichRows)),'r-'); % 2nd PC
h3=plot(squeeze(coefBoot(:,3,whichRows)),'g-'); % 3rd PC
legend([h1(1),h2(1),h3(1)],'PC1','PC2','PC3');
xlabel('Coefficient'); ylabel('Coefficient value');
title('Bootstrap replicates of 1st three PCs');

% QUESTION (Q8): Looking at the plot, what PCs are most variable? Is there a 
% pattern to this variation? 

%% Bonus: Under the hood of 'pca': Computation of the PC's.

% This section is for those wishing a more detailed explanation of what PCA
% does in terms of linear algebra operations. You can also read the
% excellent tutorial by Jon Shlens, available on the D2L course web site.

% Start with the covariance matrix:
G = cov(X);
% The PC are obtained as the eigenvectors of the covariance matrix:
[V,D] = eig(G);

% The eigenvalues give a relative measure of the variance explained by each
% PC. In MATLAB's PCA-speak, this is returned as 'latent'. But we can also
% get it as:
lambda = diag(D);   % lambda for the 1st PC is lambda(5), seemingly backwards

% For E&T p. 67-8 on interpretation of PCs:

% We can compute the first eigenvector (the first PC) by taking the last
% column of the matrix V. Transpose this into a horizontal vector. Also
% take the second eigenvector (the second to last column).

% First eigenvector:
PC1 = V(:,end)';    % [0.505,0.368,0.3456,0.451,0.535]

% TO DO (Q9): Compute the second eigenvector (PC2)
PC2 = ;

% If we wanted to best summarize each student's performance with one
% number, we would multiply each test score by the corresponding
% coefficient and sum up the result. This is called the 'loading':
loadPC1 = X(:,:) * V(:,5);  % Scientific IQ for each of the 88 students
loadPC2 = X(:,:) * V(:,4);  % loading onto 2nd PC


% QUESTION(Q10): The weights assigned by PCA can tell us about the structure of
% the multivariate data set. Interpretation of the principal components
% can sometimes be challenging. The 1st PC puts positive and roughly equal 
% weights on each of the five tests, so 'loadPC1' can be thought of as the 
% sum (or average) of each student's 5 test scores. Look at the weights in  
% PC2. What does PC2 represent? What does it mean for a student to have a 
% high'loadPC2'?


