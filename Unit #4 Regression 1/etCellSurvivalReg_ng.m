% etCellSurvivalReg.m
%
% E & T example of cell survival data, pp. 115-121
%
% RTB wrote it, Christmas eve eve, post-bikeride w/ Branfman et al.
% Converted to teaching module by SCT and RTB, August/September 2017

% Concepts covered:
% 1. quadratic regression using either 'glmfit' or 'fitglm'
% 2. unconstrained nonlinear optimization using 'fminsearch'
% 3. least median of squares regression (LMS)
% 4. regression diagnostics: Cook's distance to identify worrisome data
% 5. bootstrapping pairs to calculate SEs with LMS regression

% What to do: Login to learning catalytics and join the session for the
% module entitled "etCellSurvivalRegression, combined". You will answer a
% series of questions based on the guided programming below. Each section
% begins with a '%%'. Read through the comments and follow the instructions
% provided. In some cases you will be asked to answer a question, clearly
% indicated by 'QUESTION'. In other cases, you be asked to supply missing
% code, indicated by 'TODO'. Once you have supplied the required code, you
% can execute that section by mouse-clicking in that section (The block
% will turn yellow.) and then simultaneously hitting the 'ctrl' and 'enter'
% keys (PC) or 'command' and 'enter' keys (Mac).

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
% change the below line to redirect you to whatever this exercise is stored
% e.g. 'C://folder_name/another_folder_name/folder_where_file_is
% directory_name = ;
% cd directory_name;
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

%% Method #1: 'glmfit'

% TODO: use the function "glmfit" to fit a GLM using a model in which the
% value we retrieve is a normally distributed random variable whose mean is
% a linear function of the dose and its square (i.e. logSurvProp ~ normal
% distribution with mean = a*dose + b*dose^2, where a and b are constants).

% first make a matrix corresponding to the predictor variables (these are
% dose and dose^2)
predictors = ;

% b carries the coefficients, one for each predictor. stats carries all the
% other information about the GLM.
% HINT: use help glmfit to look up which commands to pass the function.
% Remember to set the constant term to 0!
[b,~,stats] = glmfit(predictors,ds.logSurvProp,!!! ADDITIONAL INFO HERE !!!);
b14 = b;
SE14 = stats.se;    % store the values for the std. errors of the params

% QUESTION (Q1): What do the standard errors of our parameters represent?
% Is it variability in the data, or variability in potential model fits?
% This is an important distinction!

% TODO: plot regression line using the coefficients b
ax = axis;
xVals = ax(1):ax(2);
yFit = ;
plot(xVals,yFit,'k-');

%% Method 2: 'fitglm'--gives much more information

% This method asks for more information about your data, but is
% considerably more flexible! First we define a model. Use help fitglm to
% figure out the syntax! Use the example under "formula" to figure out how
% to generate a string "modelspec" that specifies our model.

% TODO: Specify the model using a string. Again, we'll omit the intercept.
% HINT: We use -1 to tell it NOT to use a constant term. By default, it
% will use a constant term. So your modelspec string will look something
% like modelspec = logSurvProp ~ -1 + <OTHERSTUFF>; Make the names of
% variables in the string corresponds to the names of the fields in the ds
% structure!
modelspec = ;

mdl14 = fitglm(ds,modelspec,'Distribution','normal');

% mdl14 gives results in a much more user friendly format, but may be harder
% to extract info in a script. Double-click on 'mdl14' in Workspace. We see
% that mdl14 is a struct with a whole bunch of members, including,
% 'Coefficients', which is itself a table containing information on beta's,
% standard errors and t-statistics. 
beta14ls = mdl14.Coefficients.Estimate;
SE14ls = mdl14.Coefficients.SE;

% Compare the results of the two methods.
% Is the squared term justified?

% Note that we could also include an intercept in our model, either by
% including a constant term using a 1, or simply by omitting a -1 in
% modelspec.

% Also note that the 'fitglm' function also gives us much diagnostic
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

% TODO: Identify plates that may be spurious.
suspiciousPlates = ;

%% Form the fit using least median of squares (LMS)

% We will use the MATLAB function 'fminsearch':
% fminsearch finds the minimum of a scalar function of several variables,
% starting at an initial estimate. This is generally referred to as
% unconstrained nonlinear optimization.

% Type 'help optimset' to see how you control the search for a minimum
OPTIONS = optimset('Display','off','TolX',0.001);

x = ds.dose;
y = ds.logSurvProp;
% TODO: Write the function you're trying to optimize! i.e. the median of
% the squared difference between y (the logSurvProp) and your model guess
% (a constant q(1) * x plus a constant q(2) * x^2).
% HINT: Use an anonymous function so you don't have to make a separate .m
% file for this one task. fminsearch will return the values of the arguments which
% minimize the function
mdian_square_error_function = @(q) !!! FILL THIS IN;
beta14lms = fminsearch(median_square_error_function,b,OPTIONS);

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
function_to_minimize = @(q) !!! USE THE MEAN SQUARED ERROR INSTEAD OF MEDIAN;
bLS = fminsearch(function_to_minimize,bGuess,OPTIONS);

yFitLS = bLS(1).*xVals + bLS(2).*xVals.^2;
plot(xVals,yFitLS,'r--');

%% Calculate SE by bootstrapping pairs, using the LMS objective function

% Now we want to bootstrap pairs, so we'll repeat this estimation processes
nBoot = 10000;
nPts = length(ds.dose);     % number of data points in original
allBeta = zeros(nBoot,2);   % remember there are two betas

% TODO: Make this loop repeeat the estimation procedure but resampling data
% sets of this same size from our data.
for k = 1:nBoot
    bsIdx = ;   % random indexes
    x = ds.dose(bsIdx); %
    y = ds.logSurvProp(bsIdx);
    median_square_error_function = @(q) !!! Fill me in, as before.
    allBeta(k,:) = fminsearch(median_square_error_function,beta14lms,OPTIONS);
end
SE14lms = std(allBeta)';
bHat14lms = mean(allBeta);

% Question: Is our new beta estimated by repeated resampling a "better"
% estimate for the value of beta?

% NOTE: I get SE values much smaller than reported in line 3 of Table 9.5
% in E & T. But mine make more sense. When he removes the outlier (i.e.
% line 2 of table 9.5), he gets an SE for param 1 of 0.094, whereas he
% calculates the bootstrap SE for the LMS (= least median of squares) for
% param 1 to be 0.272. But why should the error be bigger for a more robust
% method? My value of 0.135 (close to the least-squares SE with the outlier removed)
% makes much more sense to me. Am I missing something? See e-mail exchange
% with Brad Efron on 12/23/2016

%% Re-do with BAD data point removed. Version #1

% TODO: Remove the suspicious data point (Plate 13) and perform the analysis
% again.
ds13 = ;  % make a copy with Plate 13 removed (try to do this through MATLAB code instead of finding the indices manually

figure, plot(ds13.dose,ds13.logSurvProp,'k+');
hold on;
xlabel('dose (rads/100)');
ylabel('log proportion alive');
title('Cell Survival Data, E&T fig. 9.3, 13 plates');

% TODO: Fit a GLM to this data using Method 1
[b,~,stats] = ;
b13 = b;
SE13 = stats.se;    % store the values for the std. errors of the params
% plot regression line
ax = axis;
xVals = ax(1):ax(2);
yFit = b(1).*xVals + b(2).*xVals.^2;
plot(xVals,yFit,'k-');

% plot the LMS fit from the full data to the data with Plate 13 removed
yFitLMS = ;
plot(xVals,yFitLMS,'k--');

%% Re-do with BAD data point removed. Version #2
figure, plot(ds.dose,ds.logSurvProp,'k+');
hold on;
plot(ds.dose(suspiciousPlates),ds.logSurvProp(suspiciousPlates),'ro');
xlabel('dose (rads/100)');
ylabel('log proportion alive');
title('Cell Survival Data, E&T fig. 9.3, 13 plates');

% Fit a GLM to the data using Method 2
mdl13 = ;
SE13ls = mdl13.Coefficients.SE;
beta13ls = mdl13.Coefficients.Estimate;

% plot regression line
ax = axis;
xVals = ax(1):ax(2);
yFit = beta13ls(1).*xVals + beta13ls(2).*xVals.^2;
h1 = plot(xVals,yFit,'r-');

yFit = beta14ls(1).*xVals + beta14ls(2).*xVals.^2;
h2 = plot(xVals,yFit,'k-');

% TODO: plot the LMS fit, too
yFitLMS = ;
h3 = plot(xVals,yFitLMS,'k--');

legend([h1,h2,h3],'Least sqares, 13 plates','Least sqares, 14 plates','Least median of squares, 14 plates');

%% Compare betas and standard errors:

% [beta14ls SE14ls beta13ls SE13ls beta14lms SE14lms]
T = table(beta14ls,SE14ls,beta13ls,SE13ls,beta14lms,SE14lms,'RowNames',{'Beta1','Beta2'})

% QUESTION: Was the plate truly spurious? What does it mean to be spurious
% in the first place? Think about which metrics can help you decide this,
% but also think about how they might sometimes fail. Arbitrary thresholds
% are always useful guides, but rarely reliable rules.