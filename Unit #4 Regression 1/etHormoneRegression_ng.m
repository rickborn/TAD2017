% etHormoneRegression.m
%
% Instructions:
% The goal of this exercise is to become familiar with the technique of
% bootstrapping and appreciate how it can be used to estimate accuracy of 
% statistics through resampling data to generate standard errors and
% confidence intervals that may otherwise be difficult to compute directly.

% What to do: Each section begins with a '%%'. Read through the comments
% and follow the instructions provided. In some cases you will be asked to
% answer a question, clearly indicated by 'QUESTION'--please do so using
% comments (i.e., by starting each line with a '%'). In other cases, you be
% asked to supply missing code, indicated by 'TODO'. Once you have supplied
% the required code, you can execute that section by mouse-clicking in that
% section (The block will turn yellow.) and then simultaneously hitting the
% 'ctrl' and 'enter' keys (PC) or 'command' and 'enter' keys (Mac).
%
% Made-up back story for neurobiologists: Prior to testing a new ALS drug
% in SOD1 mice, you need to develop and test an implantable device to allow
% chronic delivery of the drug over several days. You contract with a company
% to produce custom osmotic mini-pumps loaded with drug, and they send you
% samples of devices from three different manufacturing lots. You implant 
% them in mice for different lengths of time and then remove them and measure 
% the amount of drug remaining in the device. You would like to answer the 
% following questions:
%
% 1. How is the drug released over time?
% 2. How well might our model predict future data?
%
% Other issues to address:
% A. Different lots are likely to differ in uninteresting ways, such as the
% initial amount of drug loaded. How can we prevent this from contaminating
% our analysis?

% The data:
% Each row is data from one animal (n = 27)
% Column 1 is the manufacturing lot: A, B or C
% Column 2 is the time the device was implanted, in hours
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
ds = dataset('xlsfile','HormoneData.xlsx');

% plot it with different symbols for the different lots:
figure, gscatter(ds.hrs,ds.amount,ds.lot,'brg','xos');
hold on
%xlabel('Time worn (hrs)'); ylabel('Drug remaining (mg)');
xlabel('Time implanted (hrs)'); ylabel('Drug remaining (mg)');
%title('Drug delivery device: E & T fig. 9.1, p. 109');
title('Tests of osmotic mini-pump for drug delivery');

% QUESTION: Is the relationship between time implanted and drug remaining
% similar for the different lot numbers? By eye, does this look like a
% linear relationship?

%% Simple linear regression
% Is there a relationship between the amount of time the drug delivery
% device was implanted and the amount of drug remaining in the device? (p.
% 108). We will perform a simple linear regression.

% TO DO: Make a variable, const, which is a column of ones the length of
% one column of the dataset. We will use const to index into 
const = ones(length(ds),1);

% We will use the regress function to perform  this simple linear
% regression.
% 'regress' returns:
% estimate of beta parameters, b
% 95% CI of beta params, bCI
% residuals, r
% rejection intervals for residuals, rint (useful for dx. of outliers)
% stats: contains, in order, the R^2 statistic, the F statistic and its p-value, 
% and an estimate of the error variance
[b,bCI,r,rint,stats] = regress(ds.amount,[const,ds.hrs]);

% QUESTION: Consider the stats variable. Does the F statistic indicate a
% significant linear regression relationship between amount of time
% implanted and the amount of drug remaining in the device?

%% Simple linear regression using the GLM
mdl1 = fitglm(ds.hrs,ds.amount,'linear','Distribution','normal','Link','identity');
% NOTE: mdl1 gives results in a much more user friendly format, but may be 
% harder to extract info in a script. Double-click on 'mdl1' in Workspace.
% We see that mdl1 is a struct with a whole bunch of members, including,
% 'Coefficients', which is itself a table containing information on beta's,
% standard errors and t-statistics. 
% To get the value of beta0 (= intercept)
b0 = mdl1.Coefficients{1,'Estimate'};
b1 = mdl1.Coefficients{2,'Estimate'};

% plot the regression line
ax = axis;
xVals = [ax(1):ax(2)];
yRegression = b0 + b1.*xVals;
plot(xVals, yRegression, 'k-');
hold on;

% QUESTION: How do we answer the first question? Is the linear model a good
% predictor of drug amount? What are the values of the coefficients? (Write
% down the model!!!) What are their corresponding standard errors? 
% What is H0?

%% Linear mixed effects: allowing for different y-intercepts

% Most of the points for Lot C are above the line; most for A and B are
% below. This suggests an effect of Lot. One way to model this is with a
% linear mixed effects model. We have a fixed effect for the slope (i.e.
% 'hrs') and the intercept, but also a random addition to the intercept for
% each lot.

% E&T write this as (p. 110, eq. 9.21):
% E(y|L,x) = Beta_L + Beta_1*x
% where Beta_L takes one of three possible values (Beta_A, Beta_B, Beta_C)
% depending on which lot the data comes from.

% NOTE: MATLAB online help:
% https://www.mathworks.com/help/stats/relationship-between-formula-and-design-matrix-.html
% https://www.mathworks.com/help/stats/fitlme.html
% The formula for the model is expressed in Wilkinson notation.
% In general, a formula for model specification is a character vector of
% the form 'y ~ terms'. For the linear mixed-effects models, this formula
% is in the form 'y ~ fixed + (random1|grouping1) + ... +
% (randomR|groupingR)', where fixed and random contain the fixed-effects
% and the random-effects terms.

% TO DO:use fitlme() to make a linear mixed effect model of the amount of
% drug remaining (dependent variable) with a fixed effect of hours
% remaining and random effect of lot (see above for how to write the model
% specification. Also note that column titles are used when assigning model
% specification.)
lme = fitlme(ds,'amount ~ hrs + (1|lot)');

% Now we need to read out the individual intercepts from the model
beta = fixedEffects(lme);           % give us the fixed effects (slope & intercept)
[~,~,STATS] = randomEffects(lme);   % Compute the random-effects statistics
STATS.Level = nominal(STATS.Level); % declare as a nominal variable
betaLotA = beta(1) + STATS.Estimate(STATS.Level=='A');
betaLotB = beta(1) + STATS.Estimate(STATS.Level=='B');
betaLotC = beta(1) + STATS.Estimate(STATS.Level=='C');

% plot the individual regression lines
yRegA = betaLotA + beta(2).*xVals;
plot(xVals, yRegA, 'b-');
yRegB = betaLotB + beta(2).*xVals;
plot(xVals, yRegB, 'r-');
yRegC = betaLotC + beta(2).*xVals;
plot(xVals, yRegC, 'g-');

% QUESTION: Do the different lots have different intercepts? 
% To address the issue of lot differences, we can ask whether any of the random
% effects for the intercept are significantly different from 0. To see
% this, look at the STATS variable. What part of STATS tells you an
% intercept is significantly different from 0?

% We conclude that lot 'A' starts with significantly less drug than lots 'B' or 'C'.

%% Linear mixed effects: allowing for different intercepts AND slopes

% Make a new plot of the raw data:
figure, gscatter(ds.hrs,ds.amount,ds.lot,'brg','xos');
hold on
xlabel('Time implanted (hrs)'); ylabel('Drug remaining (mg)');
title('Tests of osmotic mini-pump for drug delivery');

% We could also allow for a random slope and intercept
lme2 = fitlme(ds,'amount ~ 1 + hrs + (1 + hrs|lot)');
beta2 = fixedEffects(lme2);           % give us the fixed effects (slope & intercept)
[~,~,STATS2] = randomEffects(lme2);   % Compute the random-effects statistics
STATS2.Level = nominal(STATS2.Level); % declare as a nominal variable

betaRand = STATS2.Estimate(STATS2.Level=='A');
yReg2A = (beta(1) + betaRand(1)) + (beta(2)+betaRand(2)).*xVals;
plot(xVals, yReg2A, 'b-');

betaRand = STATS2.Estimate(STATS2.Level=='B');
yReg2B = (beta(1) + betaRand(1)) + (beta(2)+betaRand(2)).*xVals;
plot(xVals, yReg2B, 'r-');

betaRand = STATS2.Estimate(STATS2.Level=='C');
yReg2C = (beta(1) + betaRand(1)) + (beta(2)+betaRand(2)).*xVals;
plot(xVals, yReg2C, 'g-');

% QUESTION: Does any lot elute at a greater rate than the others?
% We see that lot 'A' elutes drug at a slightly greater rate than the other 
% two. We can tell this because betaLotA is 32.3044, meaning there is less
% drug remaining after the same amount of time in lot A vs. B or C, which
% both have very similar betas.

% BONUS QUESTION: When would devices from each lot, on average, be expected
% to run out of drug? 

% Two possible approaches: 1) Write out the regression
% equation for each lot, set 'amount' to 0 and solve for x. 2) Brute force:
% xVals = 450:800;                          %Look at x values beyond the data
% yReg2C = (beta(1) + betaRand(1)) + (beta(2)+betaRand(2)).*xVals; % Redo the regression for these xVals.
% tMin = find(yReg2C== min(abs(yReg2C)));   % find y-val closest to 0
% lifeTimeC = xVals(tMin);                  % x-val corresponding to y=0

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
% empirical distribution of the residuals. 

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

% TO DO: Write the for loop that will bootstrap the residuals to allow us
% to estimate our regression coefficients (betas). We've already run 
% the regression on our data above, giving us a model (fitted data) and residuals. 
% In each iteration of the loop, we will randomly resample from the possible 
% residuals (use mdl1.Residuals.Raw to use raw residuals), giving
% us ErrorStar (see below). For each of the nBoot iterations, we will draw 
% the nPts residuals, which is the number of data points in our original dataset. 
% We will add these residuals to the fitted linear predictor of the model to 
% get yStar, which we can think of like a 'recomputation' of our data. 
% Another way to think about this is yStar(i) = BetaHat0x(i) + ErrorStar(i)
% Then, recompute the regression, only with yStar instead of the original
% Y. Store the regression coefficients from each run in a row of allBeta
% (allBeta should be a 1000x2 matrix if completed correctly).
for k = 1:nBoot
    yStar = mdl1.Fitted.LinearPredictor + mdl1.Residuals.Raw(unidrnd(nPts,nPts,1));
    [allBeta(k,:)] = regress(yStar,X);
end

% Then we can take the standard deviation of our "new" coefficients to
% estimate standard error.
bsSEresid = std(allBeta);  

% QUESTION: Compare bsSEresid with mdl1.Coefficients.SE. How similar are
% they?
% mdl1.Coefficients.SE = 0.8672, 0.0045
% bsSEresid = 0.8124, 0.0042

%% Compare with bootstrapping pairs
% Instead of resampling residuals and applying them to fitted data, we can
% select random pairs (or cases) of the data, keeping the x's and y's
% matched. We can then perform regression on these bootstrapped points.
for k = 1:nBoot
    bsIdx = unidrnd(nPts,nPts,1);
    yStar = ds.amount(bsIdx);
    xStar = [const, ds.hrs(bsIdx)];
    allBeta(k,:) = regress(yStar,xStar);
end
bsSEpairs = std(allBeta);

% QUESTION: Compare to bsSEresid and mdl1.Coefficients.SE
%% Whis is better?
% E & T state that "Bootstrapping pairs is less sensitive to assumptions
% than bootstrapping residuals." BS of residuals assumes that the error
% distribution (i.e. residuals) does not depend on x_i, and this is not
% always the case. It depends a lot on how good our assumption of linearity
% is. See fig. 9.2 on p. 114

% So regression diagnostics are important. Let's look at two measures:
% 1) Residuals vs. Fitted: 
% What we want to see: random scatter and no gross departures from linearity 
% and homoscedasticity
figure, plot(mdl1.Fitted.LinearPredictor,mdl1.Residuals.Raw,'ko');
hold on;
ax = axis;
line([ax(1),ax(2)],[0,0]);
xlabel('Linear Predictor'); ylabel('Residual');
title('Residuals vs. Fitted');

% QUESTION: What does the plot of residuals vs. fitted look like? Are our
% assumptions met?

% 2) Normal quantile plot (Q-Q Plot) of residuals
% What we want to see: points fall on main diagonal
figure, qqplot(mdl1.Residuals.Raw);

% QUESTION: What does the Q-Q Plot look like? Are our assumptions met? See
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

% Take a non-normal distribution we have good intuition about:
x = unidrnd(10,1000,1); % fat, fat tails!
figure, hist(x);
figure, qqplot(x);
% 
% % Now use qqplot to compare it to the 'right' distribution (uniform)
pd = makedist('Uniform','Lower',1,'Upper',10);
figure, qqplot(x,pd);

%% Cross-validation to measure prediction error

% Calculate the mean residual squared error of our multi-slope model
meanRSE = lme.SSE / nPts;

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

CVresiduals = zeros(nPts,1);
for k = 1:nPts
    % TO DO: specify the multi-slope model using fitlme(). Exclude the kth
    % point from the model fit. Look up fitlme() documentation for help.
    % Define the output as lmeCV.
    lmeCV = fitlme(ds,'amount ~ hrs + (1|lot)','Exclude',k);
    
    b = fixedEffects(lmeCV);           % give us the fixed effects (slope & intercept)
    [~,~,STATS] = randomEffects(lmeCV);   % Compute the random-effects statistics
    STATS.Level = nominal(STATS.Level); % declare as a nominal variable
    
    % Need to figure out which regression line to adjust slope
    switch ds.lot{k}
        case 'A'
            b(1) = b(1) + STATS.Estimate(STATS.Level=='A');
        case 'B'
            b(1) = b(1) + STATS.Estimate(STATS.Level=='B');
        case 'C'
            b(1) = b(1) + STATS.Estimate(STATS.Level=='C');
    end
    CVresiduals(k) = ds.amount(k) - (b(1) + b(2).*ds.hrs(k));
end
CVfull = sum(CVresiduals.^2) / nPts;

% QUESTION: By how many percent does meanRSE underestimate the prediction 
% error?
underEstPerCent = ((CVfull - meanRSE) / CVfull) * 100;
% 28% ! 
%% Compare residuals with CV residuals

% It looks like 'fitlme' does not automatically calculate residuals, which
% is sort of annoying. 

 lme = fitlme(ds,'amount ~ hrs + (1|lot)');
 b = fixedEffects(lme);           % give us the fixed effects (slope & intercept)
 [~,~,STATS] = randomEffects(lme);   % Compute the random-effects statistics
 STATS.Level = nominal(STATS.Level); % declare as a nominal variable
 
 realResiduals = zeros(nPts,1);
for k = 1:nPts
    % Need to figure out which regression line to use
    switch ds.lot{k}
        case 'A'
            b0 = b(1) + STATS.Estimate(STATS.Level=='A');
        case 'B'
            b0 = b(1) + STATS.Estimate(STATS.Level=='B');
        case 'C'
            b0 = b(1) + STATS.Estimate(STATS.Level=='C');
    end
    % TO DO: Write a line of code to make the kth entry of realResiduals.
    % This is the difference between the amount of drug left in the dataset
    % and the predicted amount by our model. Think about what our model
    % is here!
    realResiduals(k) = ds.amount(k) - (b0 + b(2).*ds.hrs(k));
end

% CHECK: As a reality check, if performed correctly above, myMeanRSE should
% be the same as the RSE above. 
myMeanRSE = sum(realResiduals.^2) / nPts;

figure
h1 = plot(ds.hrs,realResiduals,'ko');
hold on;
h2 = plot(ds.hrs,CVresiduals,'k*');
ax = axis;
line([ax(1),ax(2)],[0,0]);
xlabel('Time implanted (hrs)'); ylabel('residual');
title('Fig. 17.2 of E&T, p. 241: Residuals');
legend([h1,h2],'Actual residuals','Cross-validated residuals','Location','southeast');

% QUESTION: How similar are the actual and cross-validated residuals?

%% Other estimates of prediction error

% As covered on  pp. 242-243 of E&T, other measures of prediction error
% include adjusted RSE, Akaike Information Criterion (AIC) and Bayesian
% Information Criterion (BIC). We get these for free from both 'fitglm' and
% 'fitlme'. We can divide these metrics divided by the # of data points.

AIC = lme.ModelCriterion.AIC/nPts;   % Akaike Information Criterion
BIC = lme.ModelCriterion.BIC/nPts;   % Bayesian Information Criterion

% I don't get exactly same values as E&T (BIC=3.45), but they are close and
% in the right direction compared to CV. My guess is that MATLAB's are more
% accurate.

%% CV for simpler model with single slope

CVresidualsSimple = zeros(nPts,1); 
% TO DO: Consider a model with a single slope for the relationship between
% amount of drug left and time implanted, without considering differences
% between lots. Write a for loop that will specify the model, excluding the
% kth point. Then, write a line to get fixed effects of the model. Finally,
% a third line (still within the loop) will allow us to get the kth entry 
% in CVresidualsSimple by taking the difference between the amount of drug 
% left in the real data and the model's predicted value.
for k = 1:nPts
    lmeCV = fitlme(ds,'amount ~ hrs','Exclude',k);
    b = fixedEffects(lmeCV);           % give us the fixed effects (slope & intercept)
    CVresidualsSimple(k) = ds.amount(k) - (b(1) + b(2).*ds.hrs(k));
end

CVsimple = sum(CVresidualsSimple.^2) / nPts;    % E&T get 5.89; I get 6.0

% Need to get real residuals for the simple model
realResidualsSimple = zeros(nPts,1);
lmeCV=fitlme(ds,'amount ~ hrs'); 
b = fixedEffects(lmeCV); 
for k = 1:nPts
    realResidualsSimple(k)= ds.amount(k)-(b(1) + b(2).*ds.hrs(k));
end

figure
h1 = plot(ds.hrs,realResiduals,'ko');
hold on;
h2 = plot(ds.hrs,CVresiduals,'k*');
ax = axis;
h3 = plot(ds.hrs,CVresidualsSimple,'r*');
hold on
h4 = plot(ds.hrs, realResidualsSimple,'ro');
line([ax(1),ax(2)],[0,0]);
xlabel('Time implanted (hrs)'); ylabel('residual');
title('Fig. 17.2 of E&T, p. 241: Residuals');
legend([h1,h4,h2,h3],'Actual residuals (full)','Actual residuals (reduced)','CV residuals (full)',...
    'CV resiudals (reduced)','Location','southeast');

% QUESTION: How do the residuals of the data, full model, and reduced
% models compare?