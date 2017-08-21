function [ci,se] = etCorrCI(nBoot,myAlpha)

% etCorrCI.m: confidence intervals for correlation in Law School Data
%
% [ci,se] = etCorrCI(nBoot,myAlpha)
%
% Example from E & T, Table 3.2 on p. 21
%
% RTB wrote it; date unknown

% Concepts covered:
% 1. Calculating correlation coefficients with 'corr'
% 2. Bootstrapping stadard errors with built-in 'bootstrp' function
% 3. Bootstrapping confidence intervals with built-in 'bootci' function
% 4. Parametric bootstrap by sampling from a bivariate normal distribution

% cd 'C:\usr\rick\doc\Committees\PIN\PIN Director\Courses\Stats\Efron Data'

if nargin < 2, myAlpha = 0.05; end
if nargin < 1, nBoot = 5000; end

% Note: MATLAB has the random sample of 15 built in as 'lawdata'
ds82 = dataset('xlsfile','Law_School_82.xlsx'); % All law schools (census)
ds15 = dataset('xlsfile','Law_School_15.xlsx'); % random sample of 15

plot(ds82.LSAT,ds82.GPA,'k+');
hold on
plot(ds15.LSAT,ds15.GPA,'ro');
xlabel('LSAT Score'); ylabel('GPA');
% plot least squares regression line for each data set
lsline

% calculate correlation coefficients
rhoHat82 = corr(ds82.LSAT,ds82.GPA);
rhoHat15 = corr(ds15.LSAT,ds15.GPA);

% Get a bootstrap sample of correlation coefficients the old fashioned way,
% with a 'for' loop
rng default  % For reproducibility
bsRhosFL = zeros(nBoot,1);
n = length(ds15.LSAT);
for k = 1:nBoot
    thisSample = unidrnd(n,n,1);
    bsRhosFL(k) = corr(ds15.LSAT(thisSample),ds15.GPA(thisSample));
end
seBSfl = std(bsRhosFL);

% Do the same thing using built-in 'bootstrp' function
rng default  % For reproducibility
bsRhos = bootstrp(nBoot,'corr',ds15.LSAT,ds15.GPA);
figure, hist(bsRhos,30);
xlabel('Correlation coefficient'); ylabel('#');
title('Distribution of rhos from bootstrap samples of size 15')
bsAxis = axis;
seBS = std(bsRhos);
% Note: 'bootci' uses the bias correced and accelerated method by default.
% To specify method, use 
ci = bootci(nBoot,{@corr,ds15.LSAT,ds15.GPA},'alpha',myAlpha,'type','percentile');

% Now compare this to distribution obtained from random samples of size 15
% from the 'census' distribution of all 82 law schools
nTotal = length(ds82); sampSize = 15;
allSamples = unidrnd(nTotal,sampSize,nBoot);
tsRhos = diag(corr(ds82.LSAT(allSamples),ds82.GPA(allSamples)));
figure, hist(tsRhos,30);
xlabel('Correlation coefficient'); ylabel('#');
title('Distribution of rhos from true random samples of size 15 from total population of 82')
tsAxis = axis;
axis([bsAxis(1), bsAxis(2), tsAxis(3), tsAxis(4)]);

% Now compare to a 'parametric bootstrap' (p. 53 of E&T)
% "Instead of sampling with replacement from the data, we draw B samples of
% size n from the parametric estimate of the population."
% In this case, we will assume that the population has a bivariate normal
% distribution:
muHat = mean([ds15.LSAT, ds15.GPA]);
covHat = cov(ds15.LSAT,ds15.GPA);

% draw 15 pairs from this distribution and get a parametric bootstrap
% replicate of the correlation coefficient, repeat nBoot times
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

% Finally, use Fisher's transformation of rhos to get a distribution that
% should look normal with s.d. of 1/sqrt(n-3)
pbsFishersZ = 0.5 .* log((1 + pbsRhos) ./ (1 - pbsRhos));
pbsFishSE = std(pbsFishersZ)
figure, hist(pbsFishersZ,30);
xlabel('Fisher''s transformation of the correlation coefficient');
ylabel('#');

% Checking the 'untransform' calculation -- looks right
pbsRhosUT = (exp(2.*pbsFishersZ)-1) ./ (exp(2.*pbsFishersZ)+1);

% Code below is from Rebecca Senft, 7-31-2017
% BONUS: Calculate a 95% confidence interval using Fisher's transformation
% of rho. Un-transform your calculated interval using the following
% equation: untransformed = (exp(2.*Z-score)-1)/(exp(2.*Z-score)+1)
% Name the transformed upper bound 'upper' and the lower bound 'lower'. How
% does this compare to the CI we calculated when using the bootci function?
 
% Transformed CI:
boundupper = prctile(pbsFishersZ,97.5);
boundlower = prctile(pbsFishersZ,2.5);

%Converting back to untransformed units:
utupper=(exp(2.*boundupper)-1)/(exp(2.*boundupper)+1);
utlower=(exp(2.*boundlower)-1)/(exp(2.*boundlower)+1);

ciFisher=[utlower,utupper]
upperbs = prctile(bsRhos,97.5);
lowerbs = prctile(bsRhos,2.5);
ciUncorrected=[lowerbs, upperbs];

%Plot the transformed bounds on the Fisher's transformed distribution
yL = get(gca,'YLim');
line([boundupper boundupper],yL,'Color','r','LineWidth',3);
hold on
line([boundlower ;boundlower],yL,'Color','r','LineWidth',3);

% Plot the back-transformed CIs on Figure 4
figure(4);
yL = get(gca,'YLim');
line([utupper utupper],yL,'Color','r','LineWidth',2);
hold on
line([utlower;utlower],yL,'Color','r','LineWidth',2);
% compare with CIs calculated with 'bootci' and percentile method
line([ci(2) ci(2)],yL,'Color','g','LineWidth',2);
line([ci(1);ci(1)],yL,'Color','g','LineWidth',2);
% compare with CIs calculated with 'bootci' and BCA method
ciBCA = bootci(nBoot,{@corr,ds15.LSAT,ds15.GPA},'alpha',myAlpha,'type','bca');
line([ciBCA(2) ciBCA(2)],yL,'Color','y','LineWidth',2);
line([ciBCA(1);ciBCA(1)],yL,'Color','y','LineWidth',2);
