% etCellSurvivalReg.m
%
% E & T example of cell survival data, pp. 115-121
%
% RTB wrote it, Christmas eve eve, post-bikeride w/ Branfman et al.

% Concepts covered:
% 1. quadratic regression using either 'glmfit' or 'fitglm'
% 2. unconstrained nonlinear optimization using 'fminsearch'
% 3. least median of squares regression (LMS)
% 4. regression diagnostics: Cook's distance to identify worrisome data
% 5. bootstrapping pairs to calculate SEs with LMS regression

% A radiologist has run an experiment involving 14 bacterial plates. The
% plates were exposed to various doses of radiation, and the proportion of
% surviving cells measured. Larger doses lead to smaller survival
% proportions, as would be expected. There is a question mark associated
% with plate #13, as the investigator thinks it might be spurious.

% Neurobiology version, made up by RTB:
% A neurobiologist wants to test the sensitivity of a particular brain
% tumor cell line to radiation. She takes 14 plates containing tumor cells
% and exposes them to various doses of radiation, then measures the
% proportion of surviving cells. There is a question mark associated
% with plate #13, as the investigator thinks it might be spurious due to a
% malfunction of the radiation machine.

% The data:
% Each row corresponds to a plate of cells (n = 14)
% Column 1 is the plate #
% Column 2 is the radiation dose in (rads / 100)
% Column 3 is the proportion of cells surviving after radiation treatment

%% Load data and plot it

ds = dataset('xlsfile','CellSurvivalData.xlsx');
ds.logSurvProp = log(ds.survProp);  % add a column
nPts = length(ds.dose);

figure, plot(ds.dose,ds.logSurvProp,'k+');
hold on;
xlabel('dose (rads/100)');
ylabel('log proportion alive');
title('Cell Survival Data, E&T fig. 9.3, 14 plates');

%% Do least squares regression with quadratic model
% Note that we do not include an intercept term, because we know that a
% dose of zero gives survival proportion of 1, and y = log(1) = 0;

% MATLAB has a couple of different ways to fit a GLM. We'll show you both.
% One is designed to work "out of the box" whereas in the other one, you
% have much more control.

% Method #1: 'glmfit'

% TODO: use the function "glmfit" to fit a GLM using a model in which the
% value we retrieve is a normally distributed random 
[b,~,stats] = glmfit([ds.dose, ds.dose.^2],ds.logSurvProp,'normal','constant','off');
b14 = b;
SE14 = stats.se;    % store the values for the std. errors of the params
% plot regression line
ax = axis;
xVals = ax(1):ax(2);
yFit = b(1).*xVals + b(2).*xVals.^2;
plot(xVals,yFit,'k-');

% Method 2: 'fitglm'--gives much more information
modelspec = 'logSurvProp ~ -1 + dose + dose^2';
mdl14 = fitglm(ds,modelspec,'Distribution','normal')
% mdl1 gives results in a much more user friendly format, but may be harder
% to extract info in a script. Double-click on 'mdl1' in Workspace. We see
% that mdl1 is a struct with a whole bunch of members, including,
% 'Coefficients', which is itself a table containing information on beta's,
% standard errors and t-statistics. 
% To get the value of beta0 (= intercept)
% b0 = mdl14.Coefficients{1,'Estimate'};
% b1 = mdl14.Coefficients{2,'Estimate'};
beta14ls = mdl14.Coefficients.Estimate;
SE14ls = mdl14.Coefficients.SE;

% Compare the results of the two methods.
% Is the squared term justified?

% Note that we could also include an intercept in our model with:
% modelspec = 'logSurvProp ~ 1 + dose + dose^2';
% (i.e. the '-1' tells glm to omit a constant term; in fact, we don't even
% need the 1, since the default is to include a constant term. That is:
% modelspec = 'logSurvProp ~ dose + dose^2';

% Note that the 'fitglm' function also gives us much diagnostic
% information. In particular, a measure called "Cook's distance" (in
% mdl14.Diagnostics.CooksDistance) measures the effect of deleting each
% observation. A large value for a given data point, D_i, suggests that
% data point might be spurious due to either a very large residual and/or leverage.

% D_i can be interpreted as the distance one's estimates move within the
% confidence ellipsoid that represents a region of plausible values for the
% parameters. This is shown by an alternative but equivalent representation
% of Cook's distance in terms of changes to the estimates of the regression
% parameters between the cases, where the particular observation is either
% included or excluded from the regression analysis.

% Suggested cut-off values for Cook's distance are D_i > 1 or D_i > 4/n.
suspiciousPlates = ds.plateNum(mdl14.Diagnostics.CooksDistance > nPts/4);

%% Form the fit using least median of squares (LMS)

% We will use the MATLAB function 'fminsearch':
% fminsearch finds the minimum of a scalar function of several variables,
% starting at an initial estimate. This is generally referred to as
% unconstrained nonlinear optimization.

% Type 'help optimset' to see how you control the search for a minimum
OPTIONS = optimset('Display','off','TolX',0.001);

% Method #1: Old way by declaring ds as a global. Not ideal
% bLMS = fminsearch('fitFunMSR',b,OPTIONS);   % use least squares as guess

% Method #2: Using a separate function
% This is the way to do it without declaring ds as a global (thanks to Till!)
%bLMS = fminsearch(@(q)fitFunMSR2(q,ds),b,OPTIONS);

% Method #3: completely locally
x = ds.dose;
y = ds.logSurvProp;
beta14lms = fminsearch(@(q) median((y - (q(1).*x + q(2).*x.^2)).^2),b,OPTIONS);

% plot this one, too
yFitLMS = beta14lms(1).*xVals + beta14lms(2).*xVals.^2;
plot(xVals,yFitLMS,'k--');

% up-side: more robust; i.e. not pulled by outlier
% down-side: now we don't get SE for free (This is where the bootstrap
% comes in)

%% Show class that you get the exact same answer as least squares when you write 
% your own objective function. Just change 'median' to 'mean'
OPTIONS = optimset('Display','off','TolX',0.001);
x = ds.dose;
y = ds.logSurvProp;
bGuess = beta14lms;  % Try different guesses for initial betas; [0,0] works just as well
bLS = fminsearch(@(q) mean((y - (q(1).*x + q(2).*x.^2)).^2),bGuess,OPTIONS);

yFitLS = bLS(1).*xVals + bLS(2).*xVals.^2;
plot(xVals,yFitLS,'r--');

%% Calculate SE by bootstrapping pairs, using the LMS objective function
nBoot = 10000;
nPts = length(ds.dose);     % number of data points in original
allBeta = zeros(nBoot,2);   % remember there are two betas

for k = 1:nBoot
    bsIdx = unidrnd(nPts,nPts,1);   % random indexes
    x = ds.dose(bsIdx);
    y = ds.logSurvProp(bsIdx);
    allBeta(k,:) = fminsearch(@(q) median((y - (q(1).*x + q(2).*x.^2)).^2),beta14lms,OPTIONS);
end
SE14lms = std(allBeta)';
bHat14lms = mean(allBeta);  % better estimate for the betas?

% NOTE: I get SE values much smaller than reported in line 3 of Table 9.5
% in E & T. But mine make more sense. When he removes the outlier (i.e.
% line 2 of table 9.5), he gets an SE for param 1 of 0.094, whereas he
% calculates the bootstrap SE for the LMS (= least median of squares) for
% param 1 to be 0.272. But why should the error be bigger for a more robust
% method? My value of 0.135 (close to the least-squares SE with the outlier removed)
% makes much more sense to me. Am I missing something? See e-mail exchange
% with Brad Efron on 12/23/2016

%% Re-do with BAD data point removed. Version #1
ds13 = ds;  % make a copy
badID = find(ds.plateNum == 13);
if ~isempty(badID)
    ds13(badID,:) = [];
end

figure, plot(ds13.dose,ds13.logSurvProp,'k+');
hold on;
xlabel('dose (rads/100)');
ylabel('log proportion alive');
title('Cell Survival Data, E&T fig. 9.3, 13 plates');

[b,~,stats] = glmfit([ds13.dose, ds13.dose.^2],ds13.logSurvProp,'normal','constant','off');
b13 = b;
SE13 = stats.se;    % store the values for the std. errors of the params
% plot regression line
ax = axis;
xVals = ax(1):ax(2);
yFit = b(1).*xVals + b(2).*xVals.^2;
plot(xVals,yFit,'k-');

% plot the LMS fit, too
yFitLMS = beta14lms(1).*xVals + beta14lms(2).*xVals.^2;
plot(xVals,yFitLMS,'k--');

%% Re-do with BAD data point removed. Version #2
figure, plot(ds.dose,ds.logSurvProp,'k+');
hold on;
plot(ds.dose(suspiciousPlates),ds.logSurvProp(suspiciousPlates),'ro');
xlabel('dose (rads/100)');
ylabel('log proportion alive');
title('Cell Survival Data, E&T fig. 9.3, 13 plates');

modelspec = 'logSurvProp ~ -1 + dose + dose^2';
mdl13 = fitglm(ds,modelspec,'Distribution','normal','Exclude',suspiciousPlates)
SE13ls = mdl13.Coefficients.SE;
beta13ls = mdl13.Coefficients.Estimate;

% plot regression line
ax = axis;
xVals = ax(1):ax(2);
yFit = beta13ls(1).*xVals + beta13ls(2).*xVals.^2;
h1 = plot(xVals,yFit,'r-');

yFit = beta14ls(1).*xVals + beta14ls(2).*xVals.^2;
h2 = plot(xVals,yFit,'k-');

% plot the LMS fit, too
yFitLMS = beta14lms(1).*xVals + beta14lms(2).*xVals.^2;
h3 = plot(xVals,yFitLMS,'k--');

legend([h1,h2,h3],'Least sqares, 13 plates','Least sqares, 14 plates','Least median of squares, 14 plates');

%% Compare betas and standard errors:

% [beta14ls SE14ls beta13ls SE13ls beta14lms SE14lms]
T = table(beta14ls,SE14ls,beta13ls,SE13ls,beta14lms,SE14lms,'RowNames',{'Beta1','Beta2'})