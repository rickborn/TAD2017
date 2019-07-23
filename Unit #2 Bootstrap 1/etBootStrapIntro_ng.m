%% etBootStrapIntro_ng.m: Bootstrapping SE's and CI's with Grad School Data

% Instructions:
% The goal of this exercise is to become familiar with the technique of
% bootstrapping and appreciate how it can be used to estimate the precision
% of statistics through resampling data to generate standard errors and
% confidence intervals that may otherwise be difficult to compute directly.

% What to do: Login to learning catalytics and join the session for the
% module entitled "Grad School Correlations". You will answer a series of
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
% Original source of exercise:
% Efron, B. & Tibshirani Robert, J. (1993) An introduction to the
% bootstrap. Chapman & Hall, London u.a.Table 3.2 on p. 21
% 
% RTB wrote it 07 July 2019 (Kinsale, Ireland; Cork Distance Week,
% "Champion of Champions" day)
%
% see also: etCorrCI_ng_withAnswers.m

%% Concepts covered:
% 1. Standard error of the mean calculated 3 ways:
%       a) formula, b) population sampling, c) bootstrap sampling
% 2. Calculating correlation coefficients with 'corr'
% 3. Bootstrapping standard errors with the built-in 'bootstrp' function
% 4. Bootstrapping confidence intervals with the built-in 'bootci' function
% 5. Parametric bootstrap by sampling from a bivariate normal distribution
% 
% The data here are GRE (quant) and GPA (science) scores from a census of
% 82 graduate programs in neuroscience. We also have a random sample of 15
% schools from this census as well. Note that these data were collected
% prior to August of 2011, so the GRE scores were scaled from 200 to 800.

%% Define a few constants and load the data
nBoot = 10000;

% Read in the data
ds82 = readtable('Grad_School_82.xlsx'); % All graduate programs (*census*)
ds15 = readtable('Grad_School_15.xlsx'); % random sample of 15

%% Plot GPA and GRE scores
main = figure('Position',[50 10 600 900],'Name','Grad School Correlations');
subplot(3,1,1);
plot(ds82.GRE,ds82.GPA,'k+');
hold on
plot(ds15.GRE,ds15.GPA,'ro');
xlabel('GRE Score (quant)'); ylabel('GPA (science)');
% plot least squares regression line for each data set
lsline
legend('Census','Sample','Sample','Census','Location','NorthWest');

%% TODO: Calculate the mean GRE score for your sample and its standard error (SE)
% We start with something that is easy to compute directly:

meanGRE = ;
semGRE = ;

% QUESTION (Q1): What is the value of semGRE to 2 decimal places?

%% "True" standad error by sampling from the population

% TODO: Draw nBoot samples of size 15 from the CENSUS of 82, each time
% calculating the sample mean. Save each mean in 'allMeans'
nSamp = length(ds15.GRE);
nCensus = length(ds82.GRE);
allMeans = zeros(nBoot,1);

rng 'default'; % for consistency across class; You would not normally do this.
for k = 1:nBoot
    allMeans(k) = ;
end

% look at the sampling distribution of the mean
figure(1)
subplot(3,1,2);
hist(allMeans);
xlabel('mean GRE score'); ylabel('# of samples of size 15');
title('Distribution of means, sampling from census');
ax = axis;

% TODO: calculate the standard error of the mean from this sample:
semGREsamp = ;

% QUESTION (Q2): What is the value of semGREsamp to 2 decimal places?

%% Bootstrap standard error by sampling from the sample

% TODO: Calculate another SEM as you did above, but now, instead of drawing
% your samples from the CENSUS, you will draw your samples from the sample.
% You do this by sampling WITH REPLACEMENT from your original actual sample
% of the 15 graduate schools. This is the essence of the bootstrap!

allMeans = zeros(nBoot,1);
rng default;
for k = 1:nBoot
    allMeans(k) = ;
end

subplot(3,1,3);
hist(allMeans);
xlabel('mean GRE score'); ylabel('# of samples of size 15');
title('Distribution of means, re-sampling from the sample');
axis(ax);

% TODO: calculate the standard error of the mean from the bootstrap
% sampling distribution of the means:
semGREboot = ;
% QUESTION (Q3): What is the value of semGREboot to 2 decimal places?

% TODO: Compare your bootstrap estimate of the SE with that from the
% formula:
percentError = ;

% QUESTION (Q4): What is the error (in %) of the bootstrap estimate w/r/t
% that of the formula? Round to the nearest whole number in %.

%% TODO: Use the 'corr' function to calculate correlation coefficients of both the census and sample

rhoHat82 = ;
rhoHat15 = ;

% QUESTION (Q5): What is a correlation coefficient?

% QUESTION (Q6): What is the correlation coeeficient for the census? 

% QUESTION (Q7): Based on the correlation coefficient and the graph, would 
% you guess GRE score and GPA are correlated?

%% Standard error for the correlation coefficient

% Unlike for the mean, there is no handy, dandy formula for the standard
% error of a correlation coefficient. This is where the bootstrap comes in!

% Get a bootstrap sample of correlation coefficients the old fashioned way,
% using a 'for' loop

rng default  % For reproducibility
bsRhosFL = zeros(nBoot,1); 
nSamp = length(ds15.GRE);
for k = 1:nBoot
    % TO DO: Generate a variable thisSample, which will allow you to
    % randomly sample the data set by generating a list of positive whole
    % numbers of length n whose maximum value can be n. Importantly, we
    % want to sample with replacement, so thisSample also must sample with
    % replacement from 1:n). Note also that we want to preserve the
    % relationship within a given school, so we are effectively sampling
    % rows from our spreadsheet.
    thisSample = ;
    
    %Compute the correlation of GRE score ans GPA for this sample
    bsRhosFL(k) = ;
end

% Compute standard error of our correlation coefficient
seRhoBootFL = ;

% QUESTION (Q8): What is the value of seRhoBootFL to 4 decimal places?

%% Do the same thing using MATLAB's built-in 'bootstrp' function

rng default  % For reproducibility
% TO DO: Look up the documentation for 'bootstrp' and write a single line
% of code to accomplish what we did earlier with the for loop: resample
% from the sample of 15 graduate schools (ds15) nBoot times, generating a
% correlation coefficient between GRE scores and GPA each time. The output
% should be called 'bsRhos' and will be a 5000 x 1 matrix of correlation
% coefficients.
bsRhos = ;

figure
h1 = histogram(bsRhos,'Normalization','probability');
hold on;
xlabel('Correlation coefficient'); ylabel('Probability');
%title('Distribution of \rho values: re-sampling from the sample')
bsAxis = axis;
seRhoBoot = std(bsRhos);

% QUESTION (Q9): What is the mean of this distribution to 2 decimal places?
% mean(bsRhos)

% QUESTION (Q10): What is the s.d. of this distribution to 4 decimal places?

% QUESTION (Q11): How do these values compare between the 2 methods you
% used (i.e. 'for' loop vs. 'bootstrp')?

% Compare our two distributions of re-sampled correlation coefficients:
distDiff = sum(bsRhosFL - bsRhos);

% QUESTION (Q12): How can distributions using the two different methods be
% identical if they were generated by RANDOM re-sampling?

%% Sample from the census:

% As we did above for the mean, we can take advantage of the fact that we
% have data for the complete population (i.e. census), and see how our
% estimate of rho is distributed when we repeatedly sample from the
% population. That is, instead of re-sampling our sample of 15 with
% replacement, we sample the 'population' of 82 graduate schools with
% replacement.
nCensus = length(ds82.GRE);
nSamp = length(ds15.GRE);
allRhoTS = zeros(nBoot,1);

rng default
for k = 1:nBoot
    thisSample = ;
    allRhoTS(k) = ;
end

h2 = histogram(allRhoTS,'Normalization','probability');
%xlabel('Correlation coefficient'); ylabel('Probability');
%title('Distribution of \rho values: sampling from the census')
tsAxis = axis;
% axis([bsAxis(1), bsAxis(2), tsAxis(3), tsAxis(4)]);

% QUESTION (Q13): How does this distribution compare to the bootstrapped
% resampling of 15 schools? Consider the general skew, spread, location 
% (i.e. mean,median) of the distributions.

% TO DO (Q14): Compute the standard error of the correlation coefficient for
% the samples bootstrapped from the population.
seRhoBootTS = ;

%% The 'parametric bootstrap' (p. 53 of E&T)

% "Instead of sampling with replacement from the data, we draw B samples of
% size n from the parametric estimate of the population."

% The parametric bootstrap differs from the traditional bootstrap in that
% we fit a model to the data and then draw random numbers from this fitted
% model, rather than resampling the data itself. Why might one want to do
% this? Well, in rare instances when one wants to bootstrap the SE for some
% sample 'outlier', such as the 'min' or 'max', the data-driven bootstrap
% will fail. (Try this and see for yourself what is going on.) In such
% cases, the parametric bootstrap gets it right.
%
% In this case, we will assume that the population has a bivariate normal
% distribution, with means 'muHat' and a covariance matrix of 'covHat'.
muHat = mean([ds15.GRE, ds15.GPA]);
covHat = cov(ds15.GRE,ds15.GPA);

% TO DO: Using what we learned from bootstrapping by hand, create a 'for' 
% loop that uses the 'mvnrnd' function to draw nBoot samples of size sampSize 
% from a bivariate normal distribution with mean muHat and covariance covHat.
% Compute the correlation coefficient for each sample and store in a 
% variable called pbsRhos, which has been initialized for you below. 
pbsRhos = zeros(nBoot,1);
rng default
for k = 1:nBoot
    R = ;
    pbsRhos(k) = ;
end

% QUESTION (Q15): What is the standard error of the correlation coefficient
% as determined by parametric bootstrapping?
seRhoPBS = ;

h3 = histogram(pbsRhos,'Normalization','probability');
xlabel('Correlation coefficient');
% NOTE: MATLAB speaks LaTeX, so the '\rho' below will give use the Greek
% symbol for rho:
title('Distribution of \rho values')
pbsAxis = axis;
% axis([bsAxis(1), bsAxis(2), pbsAxis(3), pbsAxis(4)]);
legend([h1,h2,h3],{'Bootstrap','Census','Parametric'},'Location','Northwest');

% QUESTION (Q16): How does the SE of the correlation coefficient compare to our
% other bootstrapping strategies? If it's different, why do you think this
% may be so?

%% Side-by-side view of our histograms

figure
hist([bsRhos,allRhoTS,pbsRhos],30);
xlabel('Correlation coefficient');
ylabel('# of bootstrap replicates');
title('Distribution of \rho values')
legend('Bootstrap','Census','Parametric','Location','Northwest');
pbsAxis = axis;
axis([bsAxis(1), bsAxis(2), pbsAxis(3), pbsAxis(4)]);

%% Confidence intervals:

% We have used several different strategies to create sampling
% distributions of a given statistic:
%   1. Repeated sampling from the entire population.
%   2. Repeated re-sampling from our original sample (bootstrap)
%   3. Repeated sampling from a population defined by parameters derived
%      from our original sample (parametric bootstrap)
%
% But in each case, we have generated an estimate of the sampling
% distribution for a given statistic. Thus far, we have used these
% distributions to generate a single estimate of precision: the standard
% error. However, we can use these same distributions to calculate other
% measures of precision, such as confidence intervals. After all, under
% normal assumptions, a standard error is a kind of confidence interval,
% since we expect about 68% of the distribution to be within +/- s.d. That
% is, for example, the SEM can be thought of as defining a 68% CI for our
% estimate of the mean. But we can go further.

%% CI by asymptotic normal distribution theory

% Since the std. error is the 68% CI, we can get any other CI by just
% calculating the appropriate number of standard deviates from the normal
% distribution. Let's use our distribution of bsRhos, calculated above.
figure, 
histogram(bsRhos,'Normalization','probability');
hold on
xlabel('Correlation coefficient'); ylabel('Probability');
title('Distribution of \rho values: bootstrap')
bsAxis = axis;
axis([0.15,1.1,bsAxis(3),bsAxis(4)]);
line([mean(bsRhos),mean(bsRhos)],[bsAxis(3),bsAxis(4)],'Color','k','LineStyle','--');

% This is our mean correlation:
meanRhoBoot = mean(bsRhos);

% This is our 68% CI
seRhoBoot = std(bsRhos);

% For a 95% CI:
myAlpha = 0.05;

% You probably remember that a 95% CI is +/- 1.96 standard deviates. So we
% could calculate our CI as meanRhoBoot +/- 1.96*seRhoBoot. But say we
% wanted to be able to calculate any arbitrary confidence interval. For a
% 99% CI, we would set myAlpha to 0.01.

% TODO: Write a line of MATLAB code that will convert a desired CI,
% expressed as myAlpha to the appropriate number of standard deviates:
numStdDeviates = ;

% QUESTION (Q17): What is numStdDeviates for myAlpha = 0.001?

% TODO: Calculate the lower and upper bounds for the 95% CI
rho95CIlow = ;
rho95CIhi = ;

% QUESTION (Q18): What is the value of rho95CIhi?
% QUESTION (Q19): Does this value make sense? Why or why not?

% TODO: Draw lines for the 95%CI on our histogram in red ink:
line();
line();

%% CI by percentile method:

% In this case, we generated 10,000 samples, so a more intuitive,
% brute-force way to calculate the 95% CI is just to sort our bootstrap
% replicates and then find the values corresponding to 250th and the
% 9750th index in the sorted array.

% Sort our bootstrap replicates
bsRhosSorted = sort(bsRhos);

idxLo = ;       % index corresponding to lower bound
idxHi = ;       % index corresponding to upper bound

rho95CIpercentileLow = bsRhosSorted(idxLo);
rho95CIpercentileHi = bsRhosSorted(idxHi);

% QUESTION (Q20): What is the lower bound of the 95% CI?

% Draw lines for the 95%CI on our histogram in green ink:
line([rho95CIpercentileLow,rho95CIpercentileLow],[bsAxis(3),bsAxis(4)],...
    'Color',[0.65,0.65,0],'LineWidth',2);
line([rho95CIpercentileHi,rho95CIpercentileHi],[bsAxis(3),bsAxis(4)],...
    'Color',[0.65,0.65,0],'LineWidth',2);

%% CI using MATLAB's built-in 'bootci' function

% Note: 'bootci' uses the bias corrected and accelerated method by default.
% To specify method, indicate 'type' as 'percentile' (see help on bootci)
% TODO: Use 'bootci' to calculate the 95% CI by the percentile method:
rng default
ci = bootci();

% Draw lines for the 95%CI on our histogram in yellow ink:
% NOTE: These are identical to the ones we got from our bootstrap
line([ci(1),ci(1)],[bsAxis(3),bsAxis(4)],'Color','y');
line([ci(2),ci(2)],[bsAxis(3),bsAxis(4)],'Color','y');

% QUESTION (Q21): What is the lower bound of the 95% CI returned by
% 'bootci'?

% QUESTION (Q22): Think about this confidence interval and your earlier
% guess about whether GRE score and GPA are correlated. How can you use
% this to generate a hypothesis test? (i.e. Can we say that GRE and GPA
% are significantly correlated at p < 0.05?)

% QUESTION (Q23): Today we've explored bootstrapping as a way to estimate
% standard errors and confidence intervals for means and correlation
% coefficients. Which of these measures are the most robust across our
% different ways of bootstrapping and estimating? Which are more sensitive
% to the method we chose?