%% etCorrCI.m: Confidence intervals for correlation in Law School Data

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
% Original source of exercise:Efron, B. & Tibshirani Robert, J. (1993) 
% An introduction to the bootstrap. Chapman & Hall, London u.a.Table 3.2 on p. 21
% 
% Adapted by RTB, date unknown
% Developed for homework by RAS and RTB, July-August 2017
%% Concepts covered:
% 1. Calculating correlation coefficients with 'corr'
% 2. Bootstrapping standard errors with built-in 'bootstrp' function
% 3. Bootstrapping confidence intervals with built-in 'bootci' function
% 4. Parametric bootstrap by sampling from a bivariate normal distribution
% 
% The data here is LSAT and GPA scores from a census of 82 law schools.
% We have a random subset of 15 schools from this census as well.

%% Load data
myAlpha = 0.05;
nBoot = 5000;

% Note: MATLAB has the random sample of 15 built in as 'lawdata'
ds82 = dataset('xlsfile','Law_School_82.xlsx'); % All law schools (census)
ds15 = dataset('xlsfile','Law_School_15.xlsx'); % random sample of 15
%% Plot GPA and LSAT scores
figure
plot(ds82.LSAT,ds82.GPA,'k+');
hold on
plot(ds15.LSAT,ds15.GPA,'ro');
xlabel('LSAT Score'); ylabel('GPA');
% plot least squares regression line for each data set
lsline
legend('Census', 'Sample','Sample', 'Census','Location','NorthWest');

%% Calculate the mean LSAT for your sample and its standard error (SE)

meanLSAT = mean(ds15.LSAT);
seLSAT = std(ds15.LSAT) / sqrt(length(ds15.LSAT));

%% "True" standad error by actually sampling from the population

% draw 10000 samples of size 15 from the census of 82, each time calculating
% the sample mean
nSamp = 15;
nPop = 82;
nBoot = 10000;
allMeans = zeros(nBoot,1);

for k = 1:nBoot
    allMeans(k) = mean(ds82.LSAT(randi(nSamp,nSamp,1)));
end

figure, hist(allMeans);
seLSATsamp = std(allMeans);

% NOTE: The students might note that the standard error (SE) computed by the
% formula is a bit larger than the "true" SE calculated by sampling from
% the full population. This is a teachable moment. First, the SE we
% calculate with the formula is an *estimate*, just like the sample mean.
% We aren't guaranteed that it will be "right" any more than we are
% guaranteed that the sample mean (m) will equal the population mean (mu).
% Second, the formula for the standard error assumes that the sample size
% is much smaller than the population size, which, in this case, is not
% true. Our sample is actually about 18% of the population. In the case of
% a large sampling fraction (> 5%), we would apply a "finite population
% correction":
FPC = sqrt((nPop-nSamp)/(nPop-1));

% For more details see: Isserlis, L. (1918). "On the value of a mean as
% calculated from a sample". Journal of the Royal Statistical Society.
% 81(1): 75–81.

%% Bootstrap standard error by sampling from the sample
nSamp = 15;
nBoot = 10000;
allMeans = zeros(nBoot,1);

for k = 1:nBoot
    allMeans(k) = mean(ds15.LSAT(randi(nSamp,nSamp,1)));
end

figure, hist(allMeans);
seLSATboot = std(allMeans);

% Compare your bootstrap estimate of the SE with that from the formula

%% Compare formula SE with bootstrap SE for a bunch of different samples
nSamp = 15;
nPop = 82;
nTests = 1000;
nBoot = 5000;
allSEformula = zeros(nTests,1);
allSEboot = zeros(nTests,1);

for j = 1:nTests
    thisSample = ds82.LSAT(randi(nSamp,nSamp,1));
    allSEformula(j) = std(thisSample) / sqrt(length(thisSample));
    
    for k = 1:nBoot
        allMeans(k) = mean(thisSample(randi(nSamp,nSamp,1)));
    end
    allSEboot(j) = std(allMeans);
end

figure, plot(allSEformula,allSEboot,'b*');
hold on;
ax = axis;
line([min(ax),max(ax)],[min(ax),max(ax)]);
xlabel('SEM by formula');
ylabel('SEM by bootstrap');

meanError = mean(allSEboot ./ allSEformula);

%% Calculate correlation coefficients
rhoHat82 = corr(ds82.LSAT,ds82.GPA);
rhoHat15 = corr(ds15.LSAT,ds15.GPA);
% QUESTION (Q1): What is a correlation coefficient?

% QUESTION (Q2): How different are the correlation coefficients between the full
% census and the sample of 15 schools? 

% QUESTION (Q3): Based on the correlation coefficient and the graph, would 
% you guess LSAT score and GPA are correlated?

%% Get a bootstrap sample of correlation coefficients the old fashioned way,
% with a 'for' loop
rng default  % For reproducibility
bsRhosFL = zeros(nBoot,1); 
n = length(ds15.LSAT);
for k = 1:nBoot
    % TO DO: Generate a variable thisSample, which will allow you to
    % randomly sample the data set by generating a list of positive whole
    % numbers of length n whose maximum value can be n. Importantly, we
    % want to sample with replacement, so thisSample also must sample with
    % replacement from 1:n).
    thisSample = unidrnd(n,n,1);
    
    %Compute the correlation of LSAT score ans GPA for this sample
    bsRhosFL(k) = corr(ds15.LSAT(thisSample),ds15.GPA(thisSample));
end
% Compute standard error of our correlation coefficient
seBSfl = std(bsRhosFL);

% QUESTION (L1): Why is it that we're using the standard deviation command
% to estimate the standard error of our correlation coefficient? Why don't
% we use a standard error command here?

%% Do the same thing using built-in 'bootstrp' function
rng default  % For reproducibility

% TO DO: Look up the documentation for bootstrp and write a single line of
% code to accomplish what we did earlier with the for loop: resample from 
% the sample of 15 law schools (ds15) nBoot times, generating a correlation
% coefficient between LSAT scores and GPA each time. The output should be
% called bsRhos and will be a 5000 x 1 matrix of correlation coefficients. 
bsRhos = bootstrp(nBoot,'corr',ds15.LSAT,ds15.GPA);

figure, hist(bsRhos,30);
xlabel('Correlation coefficient'); ylabel('#');
title('Distribution of rhos from bootstrap samples of size 15')
bsAxis = axis;
seBS = std(bsRhos);

% QUESTION (Q4): Describe the distribution. What are the mean (Q3) and standard
% deviation (Q5) to four decimal points?

% We can use the bootci function to calculate a confidence interval for 
% the correlation coefficient. 
% Note: 'bootci' uses the bias corrected and accelerated method by default.
% To specify method, use 
ci = bootci(nBoot,{@corr,ds15.LSAT,ds15.GPA},'alpha',myAlpha,'type','percentile');

% QUESTION (Q6): Think about this confidence interval and your earlier
% guess about whether LSAT score and GPA are correlated. Based on CI, what
% can you say about the relationship between LSAT and GPA?
%% Now compare this to distribution obtained from random samples of size 15
% from the 'census' distribution of all 82 law schools
% Instead of re-sampling our sample of 15 with replacement, let's sample
% the 'population' of 82 law schools with replacement. 
nTotal = length(ds82); sampSize = 15;
allSamples = unidrnd(nTotal,sampSize,nBoot);

% TO DO: We need to sample correlation coefficients from the samples of 15 
% schools drawn from the census. allSamples is a 15 by 5000 matrix of
% numbers between 1 and 82 that we can use to sample the original dataset.
% Write a single line of code to compute correlation coefficients of LSAT 
% scores and GPA for these samples using allSamples to index the original
% datset. 
% HINT: You want to match indexes when computing the correlation 
% coefficients between LSAT and GPA. Corr will output a larger matrix and
% you need to extract only values from matching indices. 
tsRhos = diag(corr(ds82.LSAT(allSamples),ds82.GPA(allSamples)));

figure, hist(tsRhos,30);
xlabel('Correlation coefficient'); ylabel('#');
title('Distribution of rhos from true random samples of size 15 from total population of 82')
tsAxis = axis;
axis([bsAxis(1), bsAxis(2), tsAxis(3), tsAxis(4)]);

% QUESTION (L2): How does this distribution compare to the bootstrapped
% resampling of 15 schools? Consider the general skew, spread, location 
% (i.e. mean,median) of the distributions.

% TO DO (Q7): Compute the standard error of the correlation coefficient for
% the samples bootstrapped from the population.
seBSts =std(tsRhos)

%% Now compare to a 'parametric bootstrap' (p. 53 of E&T)
% "Instead of sampling with replacement from the data, we draw B samples of
% size n from the parametric estimate of the population."
% The parametric bootstrap differs from the traditional bootstrap in that
% we fit a model to the data and then draw random numbers from this fitted
% model, rather than resampling the data itself.
% In this case, we will assume that the population has a bivariate normal
% distribution:
muHat = mean([ds15.LSAT, ds15.GPA]);
covHat = cov(ds15.LSAT,ds15.GPA);

% TO DO: Using what we learned from bootstrapping by hand, create a for 
% loop that uses the mvnrnd function to draw nBoot samples of size sampSize 
% from a bivariate normal distribution with mean muHat andcovariance covHat.
% Compute the correlation coefficient for each sample and store in a 
% variable called pbsRhos, which has been initialized for you below. 
pbsRhos = zeros(nBoot,1);

for k = 1:nBoot
    R = mvnrnd(muHat,covHat,sampSize);
    pbsRhos(k) = corr(R(:,1),R(:,2));
end

% TO DO (Q8): What is the standard error of the correlation coefficient
% with parametric bootstrapping?
pbsSE = std(pbsRhos)

figure, hist(pbsRhos,30);
xlabel('Correlation coefficient'); ylabel('#');
title('Distribution of rhos from parametric bootstrap samples of size 15')
pbsAxis = axis;
axis([bsAxis(1), bsAxis(2), pbsAxis(3), pbsAxis(4)]);

% QUESTION (L3): How does the SE of the correlation coefficient compare to our
% other bootstrapping strategies? If it's different, why do you think this
% may be so?

% Ans: When a model fits the data properly, simulating from the model as we 
% have above generates more accurate estimates for the same sampling n than 
% re-sampling our data. Thus, the standard error for the correlation 
% coefficient is smaller here than our other bootstrapping methods. 
% HOWEVER, this assumption only works if our model is appropriate. If 
% inappropriate, we will converge on an incorrect answer. This is one 
% example of a trade-off between bias and variance. 
%% Finally, use Fisher's transformation of rhos to get a distribution that
% should look normal with s.d. of 1/sqrt(n-3)
% Fisher's transformation is useful when the population correlation 
% coefficient is not 0, which causes the sampling distribution of 
% correlation coefficients to be skewed. Without the variance-stabilizing 
% Fisher transformation, the variance would become smaller as the
% population correlation coefficient rho approaches 1. 
pbsFishersZ = 0.5 .* log((1 + pbsRhos) ./ (1 - pbsRhos));
pbsFishSE = std(pbsFishersZ)
figure, hist(pbsFishersZ,30);
xlabel('Fisher''s transformation of the correlation coefficient');
ylabel('#');
hold on

% QUESTION (Q9): Why is the standard error of the Fisher-transformed
% distribution different than our other bootstrapping distributions?

% TO DO (Q10): Compute a confidence interval for the correlation coefficient 
% with fisher's transformation of rhos. 
upper = mean(pbsRhos)+1.96.*pbsFishSE;
lower = mean(pbsRhos)-1.96.*pbsFishSE;

% QUESTION (Q11): What is the untransformed CI?
utupper=(exp(2.*upper)-1)/(exp(2.*upper)+1);
utlower=(exp(2.*lower)-1)/(exp(2.*lower)+1);
ciFish=[utlower,utupper]

% QUESTION (L4): Today we've explored bootstrapping as a means of
% estimating correlation coefficients, standard errors of rho, and CIs for
% rho. Which of these measures are the most robust across our different
% ways of bootstrapping and estimating? Which are more sensitive to the
% method we chose?