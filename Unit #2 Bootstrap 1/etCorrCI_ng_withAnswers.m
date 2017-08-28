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
% We have a subset of 15 schools from this census as well. 
%% Load data
myAlpha = 0.05
nBoot = 5000

% Note: MATLAB has the random sample of 15 built in as 'lawdata'
ds82 = dataset('xlsfile','Law_School_82.xlsx'); % All law schools (census)
ds15 = dataset('xlsfile','Law_School_15.xlsx'); % random sample of 15
%% Plot GPA and LSAT scores
plot(ds82.LSAT,ds82.GPA,'k+');
hold on
plot(ds15.LSAT,ds15.GPA,'ro');
xlabel('LSAT Score'); ylabel('GPA');
% plot least squares regression line for each data set
lsline
legend('Census', 'Sample','Sample', 'Census');

% QUESTION: Do these variables appear correlated? How strongly? How
% different is the sample from the census data?

%% Calculate correlation coefficients
rhoHat82 = corr(ds82.LSAT,ds82.GPA);
rhoHat15 = corr(ds15.LSAT,ds15.GPA);

% QUESTION: How different are the correlation coefficients between the full
% census and the sample of 15 schools?

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
% Since standard error is the standard deviation of the distribution of a 
% sample statistic, we take the standard deviation of our list of
% correlation coefficients. 
seBSfl = std(bsRhosFL);

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

% We can use the bootci function to calculate a confidence interval for 
% the correlation coefficient. 
% Note: 'bootci' uses the bias corrected and accelerated method by default.
% To specify method, use 
ci = bootci(nBoot,{@corr,ds15.LSAT,ds15.GPA},'alpha',myAlpha,'type','percentile');

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

% QUESTION: How does this distribution compare to the bootstrapped
% resampling of 15 schools? Consider the general skew, spread, location 
% (i.e. mean,median) of the distributions.

% TO DO: Compute the standard deviation of tsRhos to find the standard
% error of the correlation coefficient. 
seBSts =std(tsRhos)
% QUESTION: How does this number compare to our previous methods of 
% estimating standard error? Does this suggest our sample of 15 schools was 
% representative of the population or not?
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

pbsSE = std(pbsRhos)
figure, hist(pbsRhos,30);
xlabel('Correlation coefficient'); ylabel('#');
title('Distribution of rhos from parametric bootstrap samples of size 15')
pbsAxis = axis;
axis([bsAxis(1), bsAxis(2), pbsAxis(3), pbsAxis(4)]);

% QUESTION: What is the SE of the correlation coefficient? How does this
% compare to our other bootstrapping strategies? If it's different, why do
% you think this may be so?
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
% coefficient is greater than 0, which causes the sampling distribution of 
% correlation coefficients to be skewed. Without the variance-stabilizing 
% Fisher transformation, the variance would become smaller as the
% population correlation coefficient rho approaches 1. 
pbsFishersZ = 0.5 .* log((1 + pbsRhos) ./ (1 - pbsRhos));
pbsFishSE = std(pbsFishersZ)
figure, hist(pbsFishersZ,30);
xlabel('Fisher''s transformation of the correlation coefficient');
ylabel('#');
hold on

% QUESTION: Why is the standard error of the Fisher-transformed
% distribution different than our other bootstrapping distributions?
