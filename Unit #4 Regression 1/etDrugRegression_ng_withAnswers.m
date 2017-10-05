% etDrugRegression_ng_withAnswers.m
%
% Instructions:
% The goal of this exercise is to become familiar with the technique of
% bootstrapping and appreciate how it can be used to estimate accuracy of 
% statistics through resampling data to generate standard errors and
% confidence intervals that may otherwise be difficult to compute directly.

% What to do: Login to learning catalytics and join the session for the
% module entitled "etDrugRegression". You will answer a series of
% questions based on the guided programming below. Each section begins with
% a '%%'. Read through the comments and follow the instructions provided.
% In some cases you will be asked to answer a question, clearly indicated
% by 'QUESTION'--please do so using comments (i.e., by starting each line
% with a '%'). In other cases, you be asked to supply missing code,
% indicated by 'TODO'. Once you have supplied the required code, you can
% execute that section by mouse-clicking in that section (The block will
% turn yellow.) and then simultaneously hitting the 'ctrl' and 'enter' keys
% (PC) or 'command' and 'enter' keys (Mac).
%
% Prior to testing a new ALS drug in SOD1 mice, you need to develop and
% test an implantable device to allow chronic delivery of the drug over
% several days. You contract with a company to produce custom osmotic
% mini-pumps loaded with drug, and they send you samples of devices from
% three different manufacturing lots. You implant them in mice for
% different lengths of time and then remove them and measure the amount of
% drug remaining in the device. You would like to answer the following
% questions:
%
% 1. How is the drug released over time?
% 2. How well might our model predict future data?
%
% Other issues to address:
% Different lots are likely to differ in uninteresting ways, such as the
% initial amount of drug loaded. How can we prevent this from contaminating
% our analysis?

% The data:
% Each row is data from one animal (n = 27)
% Column 1 is the manufacturing lot: A, B or C
% Column 2 is the length of time the device was implanted, in hours
% Column 3 is the amount of drug remaining in the device, in mg.
%
% Original source of exercise: Efron, B. & Tibshirani Robert, J. (1993) 
% An introduction to the bootstrap. Chapman & Hall, London u.a. 
% Ex. 9.3 from E & T, Chapter 9, pp. 107 - 112
% Also now includes cross-validation example from Ch. 17
%
% Adapted by RTB, home with the Bubbaloo, 21 Dec. 2016
% Developed for homework by RAS and RTB, August 2017

%% Concepts covered:
% 1. plotting grouped data with 'gscatter'
% 2. simple linear regression using 'regress'
% 3. simple linear regression using 'fitglm'
% 4. linear mixed effects models using 'fitlme'
% 5. two methods for bootstrap SE estimates: residuals vs. pairs
% 6. regression diagnostics: residuals vs. fitted; q-q plot
% 7. estimating prediction error using cross-validation
%% Load and plot data
% load data, p. 107 of E&T
ds = dataset('xlsfile','DrugData.xlsx');

% plot it with different symbols for the different lots:
figure, gscatter(ds.hrs,ds.amount,ds.lot,'brg','xos');
hold on
%xlabel('Time worn (hrs)'); ylabel('Drug remaining (mg)');
xlabel('Time implanted (hrs)'); ylabel('Drug remaining (mg)');
%title('Drug delivery device: E & T fig. 9.1, p. 109');
title('Tests of osmotic mini-pump for drug delivery');

% QUESTION (Q1): By eye, does the relationship between time implanted and
% drug remaining look to be linear?

%% Simple linear regression
% Is there a relationship between the amount of time the drug delivery
% device was implanted and the amount of drug remaining in the device? (p.
% 108). We will perform a simple linear regression.

% TODO: Use the 'regress' function to perform simple linear regression.
% Look at the documentation for this function to see what you must pass to
% it and what it returns. Note that we must explicitly send a column of
% ones to represent the constant (y-intercept) for our independent
% variable:
const = ones(length(ds),1);
[betaFit,betaCI,resid,residInt,stats] = regress(ds.amount,[const,ds.hrs]);

% QUESTION (Q2): Do the confidence intervals for our beta coefficients indicate
% a significant linear regression relationship between amount of time
% implanted and the amount of drug remaining in the device?

%% Simple linear regression using the GLM

% TODO: Use the function 'fitglm' to perform the same regression:
mdl1 = fitglm(ds.hrs,ds.amount,'linear','Distribution','normal','Link','identity');

% Double-click on 'mdl1' in the Workspace. We see that mdl1 is a struct
% with a whole bunch of members, including, 'Coefficients', which is itself
% a table containing information on beta's, standard errors and
% t-statistics.

% TODO: extract the value of the y-intercept to a variable called 'b0' and
% the slope to 'b1':
b0 = mdl1.Coefficients{1,'Estimate'};
b1 = mdl1.Coefficients{2,'Estimate'};

% plot the regression line
ax = axis;
xVals = ax(1):ax(2);
yRegression = b0 + b1.*xVals;
plot(xVals, yRegression, 'k-');
hold on;

% QUESTION (Q3): What is the regression model estimate of the y-intercept? 
% (Write down the model!!!)

% QUESTION (Q4): What is its corresponding standard error?

% QUESTION (Q5): What is H0?

% QUESTION (Q6): Is our y-intercept value statistically significant at
% alpha = 0.05?

% QUESTION (Q7): Do we care?

%% Linear mixed effects: allowing for different y-intercepts

% Look closely at figure 1. Most of the points for Lot C are above the
% line; most for A and B are below. This suggests that different lots start
% out with slightly different amounts of drug. We don't really care about
% this lot-to-lot variation, but it could affect our ability to get a good
% estimate of the thing we do care about, which is the slope. To a
% statistician, variables that vary "just because" and whose individual
% labels (e.g. 'Lot A', 'Lot B', 'Lot C') don't have experimental
% significance are often referred to as "random effects," whereas variables
% that "matter" (experimentally speaking) are referred to as "fixed
% effects." Think of it this way, if the labels of the different lots got
% switched, it wouldn't much matter to us, but if the 'labels' for the
% different time points got switched, it would be a disaster. A linear
% model that allows for both 'fixed' and 'random' effects is called a
% "linear mixed effects model." We have a fixed effect for the slope (i.e.
% 'hrs') and the intercept, but we also allow a random addition to the
% intercept for each lot.

% NOTE: MATLAB online help:
% https://www.mathworks.com/help/stats/relationship-between-formula-and-design-matrix-.html
% https://www.mathworks.com/help/stats/fitlme.html
% The formula for the model is expressed in Wilkinson notation.
% In general, a formula for model specification is a character vector of
% the form 'y ~ terms'. For the linear mixed-effects models, this formula
% is in the form 'y ~ fixed + (random1|grouping1) + ... +
% (randomR|groupingR)', where fixed and random contain the fixed-effects
% and the random-effects terms.

% TO DO:use fitlme() to make a linear mixed effects model of the amount of
% drug remaining (dependent variable) with a fixed effect of hours
% remaining and random effect of lot (see above for how to write the model
% specification. Also note that column titles are used when assigning model
% specification.)
lme = fitlme(ds,'amount ~ hrs + (1|lot)');

% Now we need to read out the individual intercepts from the model
betaFit = fixedEffects(lme);           % give us the fixed effects (slope & intercept)
[~,~,STATS] = randomEffects(lme);      % Compute the random-effects statistics
STATS.Level = nominal(STATS.Level);    % declare as a nominal variable
betaLotA = betaFit(1) + STATS.Estimate(STATS.Level=='A');
betaLotB = betaFit(1) + STATS.Estimate(STATS.Level=='B');
betaLotC = betaFit(1) + STATS.Estimate(STATS.Level=='C');

% plot the individual regression lines
yRegA = betaLotA + betaFit(2).*xVals;
plot(xVals, yRegA, 'b-');
yRegB = betaLotB + betaFit(2).*xVals;
plot(xVals, yRegB, 'r-');
yRegC = betaLotC + betaFit(2).*xVals;
plot(xVals, yRegC, 'g-');

% QUESTION (Q8): Do the different lots have different intercepts? To
% address the issue of lot differences, we can ask whether any of the
% random effects for the intercept are significantly different from 0. To
% see this, look at the STATS variable. What do we conclude?

%% Compare with multiple regression using indicator variables for lot

% create indicator variables for the different lots:
lotA = strcmp(ds.lot,'A');
lotB = strcmp(ds.lot,'B');

% fit the model:
mdl1a = fitglm([ds.hrs,lotA,lotB],ds.amount,...
    'linear','Distribution','normal','Link','identity');

% When we look at the mdl1a coefficients, we see that, indeed, 'lotA' is
% significant, while 'lotB' is not.
ba0 = mdl1a.Coefficients{1,'Estimate'};
ba1 = mdl1a.Coefficients{2,'Estimate'};
ba2 = mdl1a.Coefficients{3,'Estimate'};
ba3 = mdl1a.Coefficients{4,'Estimate'};

% Re-plot the original scatter plot with the simple regression line:
figure, gscatter(ds.hrs,ds.amount,ds.lot,'brg','xos');
hold on
xlabel('Time implanted (hrs)'); ylabel('Drug remaining (mg)');
title('Tests of osmotic mini-pump for drug delivery');
ax = axis;
xVals = ax(1):ax(2);
yRegression = b0 + b1.*xVals;
plot(xVals, yRegression, 'k-');

% Plot the regression line for lotA:
yVals = ba0 + (ba1 .* xVals) + (ba2 .* ones(size(xVals)));
plot(xVals,yVals,'b-');

% Plot the regression line for lotB:
yVals = ba0 + (ba1 .* xVals) + (ba3 .* ones(size(xVals)));
plot(xVals,yVals,'r-');

% Plot the regression line for lotC:
yVals = ba0 + (ba1 .* xVals);
plot(xVals,yVals,'g-');

% We get the same answers as we did with the LME model!

%% Application of the bootstrap (p. 111): bootstrap the residuals

% Classic quote: "Thus reassured that the bootstrap is giving reasonable
% answers in a case we can analyze mathematically, we can go on to apply
% the bootstrap to more general regression models that have no mathematical
% solution: where the regression function is non-linear in the parameters
% beta, and where we use fitting methods other than least-squares."

% Basic idea is that we need an estimate of both the regression
% coefficients (beta) and the PDF of the error terms (F). So we use our
% estimate of beta to calculate the approximate errors, e_i = y_i - BX.
% These are just the residuals. In other words, our estimate of F is the
% empirical distribution of the residuals!

% Recall that from our simple regression with the GLM, we get estimates of the
% standard error for each of our coefficients:
mdl1.Coefficients.SE;
% intercept SE = 0.8672
% slope SE = 0.0045

% We also get everything we need for the bootstrap:
% The fitted values (i.e. the predicted values for our actual x-values):
figure(1)
hP = plot(ds.hrs,mdl1.Fitted.LinearPredictor,'ko');
set(hP,'MarkerSize',3,'MarkerFaceColor','k');

% . . . and the residuals (in 4 different flavors; 
nBoot = 1000;
nPts = length(ds.hrs);  % number of data points in original
allBeta = zeros(nBoot,2);
X = [const,ds.hrs];

% TO DO: Write the 'for' loop that will bootstrap the residuals to allow us
% to estimate our regression coefficients (betas). We've already run 
% the regression on our data above, giving us a model (fitted data) and residuals. 
% In each iteration of the loop, we will randomly resample from the possible 
% residuals (use mdl1.Residuals.Raw to use raw residuals), giving
% us errorStar (see below). For each of the nBoot iterations, we will draw 
% the nPts residuals, which is the number of data points in our original dataset. 
% We will add these residuals to the fitted linear predictor of the model to 
% get yStar, which we can think of like a 'recomputation' of our data. 
% Another way to think about this is yStar(i) = betaHat0x(i) + errorStar(i)
% Then, recompute the regression, only with yStar instead of the original
% Y. Store the regression coefficients from each run in a row of allBeta
% (allBeta should be a nBoot-x-2 matrix if completed correctly). We want to
% use 'regress' within our 'for' loop, because 'fitglm' is much slower (It
% is calculating a bunch of stuff that we don't need in the bootstrap.)
rng default
for k = 1:nBoot
    yStar = mdl1.Fitted.LinearPredictor + mdl1.Residuals.Raw(unidrnd(nPts,nPts,1));
    [allBeta(k,:)] = regress(yStar,X);
end

% Then we can take the standard deviation of our "new" coefficients to
% estimate standard error.
bsSEresid = std(allBeta);  

% THOUGHT QUESTION (no LC component): Compare bsSEresid with
% mdl1.Coefficients.SE. How similar are they?
% mdl1.Coefficients.SE = 0.8672, 0.0045
% bsSEresid = 0.8542, 0.0043

% QUESTION (Q9): What is your estimate of the standard error of the slope
% coefficient based on bootstrapping residuals?

%% Bootstrapping with pairs

% Instead of resampling residuals and applying them to fitted data, we can
% select random pairs (or cases) of the data, keeping the x's and y's
% matched. That is, if we picture the row (observations) by columns
% (variables) structure of our data, we are randomly sampling rows.
% We then perform regression on these bootstrapped pairs. Again, use the
% 'regress' function.

rng default
for k = 1:nBoot
    bsIdx = unidrnd(nPts,nPts,1);
    yStar = ds.amount(bsIdx);
    xStar = [const, ds.hrs(bsIdx)];
    allBeta(k,:) = regress(yStar,xStar);
end
bsSEpairs = std(allBeta);

% TODO: Compare bsSEpairs to bsSEresid and mdl1.Coefficients.SE
% mdl1.Coefficients.SE = 0.8672, 0.0045
% bsSEresid = 0.8542, 0.0043
% bsSEpairs = 0.7781, 0.0043

% QUESTION (Q10): What is your estimate of the standard error of the slope
% coefficient based on bootstrapping pairs?

%% Which is better?
% E & T state that "Bootstrapping pairs is less sensitive to assumptions
% than bootstrapping residuals." BS of residuals assumes that the error
% distribution (i.e. residuals) does not depend on x_i, and this is not
% always the case. It depends a lot on how good our assumption of linearity
% is and on the homoscedasticity of the data. See fig. 9.2 on p. 114 of E&T

% So regression diagnostics are important. Let's look at two measures:
% 1) Residuals vs. Fitted: 
% What we want to see: random scatter and no gross departures from linearity 
% and homoscedasticity.
figure, plot(mdl1.Fitted.LinearPredictor,mdl1.Residuals.Raw,'ko');
hold on;
ax = axis;
line([ax(1),ax(2)],[0,0]);
xlabel('Linear Predictor'); ylabel('Residual');
title('Residuals vs. Fitted');

% THOUGHT QUESTIONS: What does the plot of residuals vs. fitted look like? Are our
% assumptions met?

%QUESTION (Q11): What do we mean by 'homoscedasticity'?

% 2) Normal quantile plot (Q-Q Plot) of residuals
% What we want to see: points fall on main diagonal
figure, qqplot(mdl1.Residuals.Raw);

% THOUGHT QUESTIONS: What does the Q-Q Plot look like? Are our assumptions met? See
% the next section for a better intuition on what Q-Q plots look like with
% matching vs. non-matching distributions.

%% Bonus on intuition for Q-Q Plot
% from MATLAB' documentation:
% A quanitle-quantile plot (also called a q-q plot) visually assesses
% whether sample data comes from a specified distribution. Alternatively, a
% q-q plot assesses whether two sets of sample data come from the same
% distribution.
% 
% A q-q plot orders the sample data values from smallest to largest, then
% plots these values against the expected value for the specified
% distribution at each quantile in the sample data. The quantile values of
% the input sample appear along the y-axis, and the theoretical values of
% the specified distribution at the same quantiles appear along the x-axis.
% If the resulting plot is linear, then the sample data likely comes from
% the specified distribution.
% 
% The q-q plot selects quantiles based on the number of values in the
% sample data. If the sample data contains n values, then the plot uses n +
% 1 quantiles. Plot the ith ordered value (also called the ith order
% statistic) against the i (n+1) th quantile of the specified distribution.
% 
% A q-q plot can also assesses whether two sets of sample data have the
% same distribution, even if you do not know the underlying distribution.
% The quantile values for the first data set appear on the x-axis and the
% corresponding quantile values for the second data set appear on the
% y-axis. Since q-q plots rely on quantiles, the number of data points in
% the two samples does not need to be equal. If the sample sizes are
% unequal, the q-q plot chooses the quantiles based on the smaller data
% set. If the resulting plot is linear, then the two sets of sample data
% likely come from the same distribution.

% Start with a non-normal distribution we have good intuition about:
% TODO: Take 10,000 draws from a uniform discrete random distribution with
% a maximum of 10 and store it in a variable called 'A'
A = unidrnd(10,10000,1);
figure, subplot(3,1,1)
histogram(A);

% TODO: Make a q-q plot vs. the percentiles in a normal distribution
subplot(3,1,2);
qqplot(A);
 
% TODO: Now use qqplot to compare it to the 'right' distribution (uniform)
pd = makedist('Uniform','Lower',1,'Upper',10);
subplot(3,1,3);
qqplot(A,pd);

%% Cross-validation to measure prediction error

% Calculate the mean residual squared error of our original, simple model:
meanRSE = sum(resid.^2) / nPts;

% The meanRSE is a measure of how well our model describes the data. But it
% is overly optimistic, because it is measuring performance using the same
% data that was used to fit the model. That is, the model is optimized to
% fit precisely this data. But how well would it do on another data set?
% Well, we could either repeat the experiment, and see how well the model
% fit to the first data set predicted the new experimental values.
% Alternatively, we can use cross-validation, in which we divide up the
% data into K equal-sized parts, fit the model to the other K-1 parts, then
% measure the error in predicting the Kth part. This is called K-fold
% cross-validation. A common method is K=n, referred to as "leave-one-out"
% cross-validation.

% TO DO: Perform a leave-one-out cross-validation in order to calculate our
% 'CVresiduals':
CVresiduals = zeros(nPts,1);
for k = 1:nPts
    sV = true(nPts,1);      % define a selection vector
    sV(k) = 0;
    [betaFit] = regress(ds.amount(sV),[const(sV),ds.hrs(sV)]);
    CVresiduals(k) = ds.amount(k) - (betaFit(1) + betaFit(2).*ds.hrs(k));
end

% TODO: Calculate the mean squared error of our cross-validated residuals.
CVfull = sum(CVresiduals.^2) / nPts;

% QUESTION (Q12): By how many percent does meanRSE underestimate the prediction 
% error as computed by cross-validation? Use 'CVfull' as the gold standard.
underEstPerCent = ((CVfull - meanRSE) / CVfull) * 100;  % 13%

%% Compare actual residuals with the residuals obtained by cross-validation

figure
h1 = plot(ds.hrs,resid,'ko');
hold on;
h2 = plot(ds.hrs,CVresiduals,'k*');
ax = axis;
line([ax(1),ax(2)],[0,0]);
xlabel('Time implanted (hrs)'); ylabel('residual');
title('Real vs. CV Residuals');
legend([h1,h2],'Actual residuals','Cross-validated residuals','Location','southeast');

% THOUGHT QUESTION: How similar are the actual and cross-validated residuals?

%% Other estimates of prediction error

% As covered on  pp. 242-243 of E&T, other measures of prediction error
% include adjusted RSE, Akaike Information Criterion (AIC) and Bayesian
% Information Criterion (BIC). 

% TODO: Read up on 'Akaike Information Criterion', then figure out how to
% derive these values for our two models. HINT: We get these for free from
% both 'fitglm' and 'fitlme'.

AICsimple = mdl1.ModelCriterion.AIC;
BICsimple = mdl1.ModelCriterion.BIC;
AIClme = lme.ModelCriterion.AIC;   
BIClme = lme.ModelCriterion.BIC;

% QUESTION (Q13): Based on the AIC values, which is the "better" model?
