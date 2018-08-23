% etCellSurvivalReg_ng.m
%
% E & T example of cell survival data, pp. 115-121
%
% RTB wrote it, Christmas eve eve, post-bikeride w/ Branfman et al.
% Converted to teaching module by SCT and RTB, August/September 2017

% Concepts covered:
% 1. quadratic regression using either 'glmfit' or 'fitglm'
% 2. unconstrained nonlinear optimization using 'fminsearch'
% 3. anonymous functions
% 4. least median of squares regression (LMS)
% 5. regression diagnostics: Cook's distance to identify worrisome data
% 6. bootstrapping pairs to calculate SEs with LMS regression

% What to do: Login to Learning Catalytics (LC) and join the session for
% the module entitled "Cell Survival Regression". You will answer
% a series of questions based on the guided programming below. Each section
% begins with a '%%'. Read through the comments and follow the instructions
% provided. In some cases you will be asked to answer a question, clearly
% indicated by 'QUESTION' and a corresponding 'Q#' that directs you to
% answer the relevant question in LC. In other cases, you be asked to
% supply missing code, indicated by 'TODO'. Once you have supplied the
% required code, you can execute that section by mouse-clicking in that
% section (The block will turn yellow.) and then simultaneously hitting the
% 'ctrl' and 'enter' keys (PC) or 'command' and 'enter' keys (Mac).

% A neurobiologist wants to test the sensitivity of a particular brain
% tumor cell line to radiation. She takes 14 plates containing tumor cells
% and exposes them to various doses of radiation, then measures the
% proportion of surviving cells. There is a question mark associated
% with plate #13, as the investigator thinks it might be spurious due to a
% malfunction of the radiation machine. This was clearly written in her lab
% notebook at the time of the experiment.

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
title('Cell Survival Data, 14 plates');

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

% First, make a matrix corresponding to the predictor variables (these are
% dose and dose^2)
predictors = ;

% b14 carries the coefficients, one for each predictor. stats carries all
% the other information about the GLM. HINT: use 'doc glmfit' to look up
% which commands to pass the function. Remember to set the constant term to
% 0.
[b14,~,stats] = glmfit(predictors,ds.logSurvProp,!!! ADDITIONAL INFO HERE !!!);
SE14 = stats.se;    % store the values for the std. errors of the params

% TODO: plot regression line using the coefficients b14
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
% structure.
modelspec = ' ';

mdl14 = fitglm(ds,modelspec,'Distribution','normal');

% mdl14 gives results in a much more user friendly format, but may be harder
% to extract info in a script. Double-click on 'mdl14' in Workspace. We see
% that mdl14 is a struct with a whole bunch of members, including,
% 'Coefficients', which is itself a table containing information on beta's,
% standard errors and t-statistics. 
beta14ls = mdl14.Coefficients.Estimate;
SE14ls = mdl14.Coefficients.SE;

% Compare the results of the two methods. They should obviously give you
% the same answer for the same model with the same data.

% QUESTION (Q1): Is the squared term justified?

% QUESTION (Q2): What average log survival proportion does your
% least-squares regression GLM fit make for a radiation dose of 6 rads /
% 100?

% QUESTION (Q3): What is the standard error for the quadratic term?

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

% Recall that standard linear regression minimizes the squared distance of
% the y-value of each data point from the regression line. Using this
% metric has many mathematical advantages, but there are other so-called
% "loss functions" that we might choose. One disadvantage of least-squares
% is that it is very sensitive to outliers: its distance is blown up by the
% squaring, and then it disproportionately influences the mean. But we know
% a more robust statistic in such circumstances: the median. Since it is
% the mid-point of the sorted data, one large outlier won't have much
% influence on the median. So what we'll do here is to essentially write
% our own loss function, then use a method called unconstrained nonlinear
% optimization to find the model parameters that minimize this function.
% The MATLAB function that performs this optimization is called
% 'fminsearch'.

% Type 'help optimset' to see how you control the search for a minimum
OPTIONS = optimset('Display','off','TolX',0.001);

% This just makes it easier to fit our function on one line:
x = ds.dose;
y = ds.logSurvProp;

% TODO: Write the function you want to optimize. i.e. the median of the
% squared difference between y (the logSurvProp) and your model guess (a
% constant q(1) * x plus a constant q(2) * x^2).
% HINT: Use an anonymous function so you don't have to make a separate .m
% file for this one task. fminsearch will return the values of the
% arguments which minimize the function
median_square_error_function = @(q) !!! FILL THIS IN;

% NOTE: 'fminsearch' requires that we supply an initial guess as to the
% values of our parameters (this is the 2nd argument passed to
% 'fminsearch'). Since there is no guarantee that 'fminsearch' will find a
% global minimum, it may be important to start with a decent guess. Since
% we have already done least-squares regression, we will use these values
% as our initial guess for the parameters ('b14').
beta14lms = fminsearch(median_square_error_function,b14,OPTIONS);

% plot this one, too
yFitLMS = beta14lms(1).*xVals + beta14lms(2).*xVals.^2;
plot(xVals,yFitLMS,'k--');

% QUESTION (Q4): Does using the least median square error significantly
% change the values of the coefficients?

% QUESTION (Q5): Why is your least median square error fit so different
% from least mean square error?

%% Use 'fminsearch' to perform standard least squares regression

% Instead of using 'glmfit', we should be able to achieve the same result
% by using 'fminsearch' with the appropriate function to minimize (referred
% to as an "objective function").

% QUESTION (Q6): What do we need to change in our above objective function
% ('median_square_error_function') to perform standard least-squares
% regression?
OPTIONS = optimset('Display','off','TolX',0.001);
x = ds.dose;
y = ds.logSurvProp;
bGuess = beta14lms;  % Try different guesses for initial betas;

% TODO: Write an objective function for least-squares regression
function_to_minimize = @(q) !!! USE THE MEAN SQUARED ERROR INSTEAD OF MEDIAN;

% Perform the minimization to find optimal values of the model parameters:
bLS = fminsearch(function_to_minimize,bGuess,OPTIONS);

yFitLS = bLS(1).*xVals + bLS(2).*xVals.^2;
plot(xVals,yFitLS,'r--');

% QUESTION (Q7): Compare the beta coefficients you got with this procedure
% with those from 'glmfit'. Are they identical?

% QUESTION (Q8): If the answer to the above question is 'no', why are they
% not identical?

%% Calculate SE by bootstrapping pairs, using the LMS objective function

% Now we want to bootstrap pairs, so we'll repeat this estimation processes
nBoot = 10000;
nPts = length(ds.dose);     % number of data points in original
allBeta = zeros(nBoot,2);   % remember there are two betas

% TODO: Make this loop repeeat the estimation procedure but resampling data
% sets of this same size from our data.
rng default
for k = 1:nBoot
    bsIdx = ;   % random indexes
    x = ds.dose(bsIdx); %
    y = ds.logSurvProp(bsIdx);
    median_square_error_function = @(q) !!! Fill me in, as before.
    allBeta(k,:) = fminsearch(median_square_error_function,beta14lms,OPTIONS);
end
SE14lms = std(allBeta)';    % transpose so that it has the same form as other SEs
bHat14lms = mean(allBeta);

% QUESTION (Q9): What is your estimate for the quadratic term from taking
% the mean across your bootstrapped samples?

% QUESTION (Q10): Is our new beta estimated by repeated resampling a
% "better" estimate for the value of beta?

% TODO: Calculate the 95% CI for the quadratic term using the percentile
% method:
!!! Your code here

% QUESTION (Q11): Based on the bootstrap, does the quadratic term belong in
% the model?

%% Re-do standard least-squares with suspicious data point removed.

% You may recall that when we used Method #2 ('fitglm') to do our
% regression, we used a measure called "Cook's distance" to identify plates
% that were suspicious on statistical grounds. You might also recall that
% we were suspicious, a priori, about plate 13 due to a comment in the
% experimnter's notebook. So let's see what our standard least-squares
% regression looks like if we omit this data point.

% TODO: Remove the suspicious data point (Plate 13) and perform the
% analysis again. Make a copy of the data with Plate 13 removed (try to do
% this through MATLAB code instead of finding the indices manually
ds13 = ;  

% Make a new figure and plot all of the data, but with the suspicious plate
% also circled in red ink:
figure, plot(ds.dose,ds.logSurvProp,'k+');
hold on;
plot(ds.dose(suspiciousPlates),ds.logSurvProp(suspiciousPlates),'ro');
xlabel('dose (rads/100)');
ylabel('log proportion alive');
title('Cell Survival Data');

% Fit our basic quadratic model to ds13 using 'fitglm'
mdl13 = ;
SE13ls = mdl13.Coefficients.SE;
beta13ls = mdl13.Coefficients.Estimate;

% Plot the regression line fit to ds13
ax = axis;
xVals = ax(1):ax(2);
yFit = beta13ls(1).*xVals + beta13ls(2).*xVals.^2;
h1 = plot(xVals,yFit,'r-');

% Re-plot our original regression line to the full data set
yFit = beta14ls(1).*xVals + beta14ls(2).*xVals.^2;
h2 = plot(xVals,yFit,'k-');

% TODO: Plot the fit made using least MEDIAN squares to the full data set.
yFitLMS = ;
h3 = plot(xVals,yFitLMS,'k--');

legend([h1,h2,h3],'Least sqares, 13 plates','Least squares, 14 plates','Least median of squares, 14 plates');

%% Compare betas and standard errors:

% This just formats our results into something easier to grok:
T = table(beta14ls,SE14ls,beta13ls,SE13ls,beta14lms,SE14lms,'RowNames',{'Beta1','Beta2'})

% THOUGHT QUESTIONS (No LC component): Was plate #13 truly spurious? What
% does it mean to be spurious in the first place? Think about which metrics
% can help you decide this, but also think about how they might sometimes
% fail. Arbitrary thresholds are always useful guides, but rarely reliable
% rules.