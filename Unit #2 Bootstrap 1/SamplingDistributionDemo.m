% SamplingDistributionDemo.m

% Start with the left half of E&T's diagram: "Real World"
% Why do many of the things we measure have approximately normal
% distributions? (e.g. hemoglobin conc. of Swedish women, heights of boys
% in the 5th grade, measurement errors)
%
% Ans. Many natural phenomena are the result of a bunch of random events,
% like coin tosses, and the CLT tells us that such sums (or means) are
% normally distributed.

%% Distribution of binomial sums: CLT
% Flip a fair coin n times and count the number of heads

Ns = [2,3,5,7,10,20];
h1 = figure('Name','Flipping coins and counting heads');
for k = 1: length(Ns)
    A = sum(round(rand(Ns(k),10000)));
    figure(h1);
    subplot(3,2,k);
    hist(A,[0:Ns(k)]);
    xlabel('# of Heads');
    ylabel('#');
    title([num2str(Ns(k)) ' Coin Tosses']);
end

%% Simple version of sums from uniform discrete random distribution
% Example: rolling a die

% What is the expectation?
muDie = sum([1:6] .* (1/6));    % 3.5, (n+1)/2

% What is the variance?
varDie = sum(([1:6].^2) .* (1/6)) - (3.5).^2;   % 2.9167, (n^2 - 1)/12

% Let's roll a die 10 times and take the sum
nRolls = 10;
nSims = 10000;
allRolls = unidrnd(6,nRolls,nSims);

figure, histogram(sum(allRolls),'Normalization','pdf');
hold on;
xlabel('sum of 10 rolls of a die');
ylabel('probability');

% How well does this match the normal approximation?
ax = axis;
xVals = ax(1):0.1:ax(2);
yVals = normpdf(xVals,muDie*nRolls,sqrt(varDie*nRolls));
plot(xVals,yVals,'r-','LineWidth',2);

% Look at this using qqplot
figure, qqplot(sum(allRolls));

%% Similar but with coin tosses
% What is the expectation?
% Let's flip a coin 10 times and take the sum of the # of heads
nFlips = 10;
nSims = 10000;
pHead = 0.5;
allFlips = round(rand(nFlips,nSims));

% What is the expectation?
muCoin = nFlips * pHead;

% What is the variance?
varCoin = nFlips * pHead * (1 - pHead);

figure, histogram(sum(allFlips),'Normalization','pdf');
hold on;
xlabel('sum of # of heads on 10 flips of a fair coin');
ylabel('probability');

% How well does this match the normal approximation?
ax = axis;
xVals = ax(1):0.1:ax(2);
yVals = normpdf(xVals,muCoin,sqrt(varCoin));
plot(xVals,yVals,'r-','LineWidth',2);

% Look at this using qqplot
figure, qqplot(sum(allFlips));

%% When is the normal a reasonable approximation to the binomial?
Ns = [2,4,8,10:10:1000];
p = 0.5;

binoProbs = zeros(size(Ns));
normProbs = zeros(size(Ns));

% calculate probability of obtaining the mean (n*p)
% . . . or one of the extremes (e.g. n)
for k = 1:length(Ns)
    n = Ns(k);
    binoProbs(k) = binopdf(n*p,n,p);
    normProbs(k) = normpdf(n*p,n*p,sqrt(n*p*(1-p)));
end
figure, plot(binoProbs,normProbs,'k.');
hold on;
ax = axis;
h1 = line([0,max(ax)],[0,max(ax)]);
set(h1,'Color','r','LineStyle','--');
xlabel('Binomial Probability')
ylabel('Normal Probability');

figure, loglog(Ns,normProbs ./ binoProbs,'k*');
xlabel('N');
ylabel('normpdf / binopdf');

%% How do we tell if something is normally distributed?

% Visual impression: Q-Q plot
figure, qqplot(A);

% see also: probplot

% Statistical tests:
% Lilliefors Test: lillietest
% Jarque-Bera test: jbtest
% One-sample Kolmogorov-Smirnov test: kstest
% Anderson-Darling test: adtest
% Shapiro-Wilk goodness-of-fit test for normality: download from MathWorks
%
% The rap against most of these tests is that they are not very powerful,
% which, in this case, means that they are not very conservative as
% deciders of normality. How would you get a sense of this by simulation?


%% What is the power and beauty of 'knowing' something is normally distributed?
%  heights of 99 five-year-old British boys in cm
load boysHgtsCm

%figure('Name','Heights of 99 five-year-old British boys in cm');
% subplot(2,2,1)
% hist(H,2);
% xlabel('Height (cm)');
% ylabel('#');
% subplot(2,2,2)
% hist(H,6);
% xlabel('Height (cm)');
% ylabel('#');
% subplot(2,2,3)
% hist(H,15);
% xlabel('Height (cm)');
% ylabel('#');
% subplot(2,2,4)
% hist(H,50);
hf = figure;
h = histogram(H,'Normalization','pdf');
% Probability density function estimate. The height of each bar is, (number
% of observations in the bin) / (total number of observations * width of
% bin). The area of each bar is the relative number of observations. The
% sum of the bar areas is 1.
hold on;
xlabel('Height (cm)');
ylabel('Probability');
title('Heights of 99 five-year-old British boys in cm');

% To normalize on our own:
figure;
h = histogram(H);
pdf = h.Values ./ (h.BinWidth * sum(h.Values));
binCtrs = h.BinEdges(1:end-1) + h.BinWidth/2;
figure, bar(binCtrs,pdf);

% Jarque-Bera test for normality; 1 means we reject H0 that data is
% normally distributed
figure, qqplot(H);
%[h,p] = jbtest(H);

% Compare with t-distribution
% nBoys = length(H);
% B = 10000;
% ps = linspace(1/(B+1),1-1/(B+1),B);
% tPercentiles = tinv(ps,length(H)-1);
% figure, qqplot(H,tPercentiles);

% The true beauty and power is that, now, knowing just these two numbers, I
% can tell you the probability that a given boy is in soome height range.
% E.g. What is the probability that a boy is between 100 and 110 cm?
muHgt = mean(H);
sigmaHgt = std(H);
pHgt = normcdf(110,muHgt,sigmaHgt) - normcdf(100,muHgt,sigmaHgt);

% superimpose nl distribution on histogram
figure(hf);
ax = axis;
xVals = ax(1):0.01:ax(2);
yVals = normpdf(xVals,muHgt,sigmaHgt);
plot(xVals,yVals,'r-');

% calculate the SEM in two ways: formula and bootstrap
n = length(H);
semHgt = std(H) / sqrt(n);
% Ans = 0.5232

nBoot = 10000;
allMuHgt = zeros(nBoot,1);
n = length(H);
rng default
for k = 1:nBoot
    allMuHgt(k) = mean(H(unidrnd(n,n,1)));
end
semHgtBS = std(allMuHgt);
% Ans = 0.5216

% or use MATLAB's 'bootstrp' function
rng default
[bootStat] = bootstrp(nBoot,@mean,H);
std(bootStat)
% Ans = 0.5216

%% Distribution of the t-statistic
Ns = [5:5:30];
B = 10000;
ps = linspace(1/(B+1),1-1/(B+1),B);
h1 = figure('Name','Comparison with the t distribution');
h2 = figure('Name','Comparison with the normal distribution');

for k = 1:length(Ns)
    
    % draw n samples from N(0,1)
    X = randn(Ns(k),B);
    
    % Can also compare what happens when we use the sample estimate of the
    % standard deviation vs. the known SD of 1. With the former, the
    % t-distribution is good and normal is bad; vice versa when we
    % normalize with known SD of 1.
    % t-statistic is the mean normalized by the standard error
    tStats = sqrt(Ns(k)) .* (mean(X) ./ std(X));
    %tStats = sqrt(Ns(k)) .* (mean(X) ./ 1);
    
    % compare with t distribution
    figure(h1);
    subplot(3,2,k);
    tPercentiles = tinv(ps,Ns(k)-1);
    qqplot(tStats,tPercentiles);
    xlim([-4.5,4.5]);
    title(['N = ' num2str(Ns(k))]);
    
    % compare with standard normal distribution
    figure(h2)
    subplot(3,2,k);
    qqplot(tStats)
    xlim([-4.5,4.5]);
end

%% Compare sampling distributions for the mean and the median
Ns = [5:10:55];
B = 10000;
meanSD = zeros(size(Ns));
medianSD = meanSD;
figure
for k = 1:length(Ns)
    X = randn(Ns(k),B);
    allMu = mean(X);
    allMed = median(X);
    subplot(3,2,k)
    hist([allMu',allMed'],25)
    legend('Means','Medians');
    title(['N = ' num2str(Ns(k))]);
    meanSD(k) = std(allMu);
    medianSD(k) = std(allMed);
end
figure
plot(Ns,meanSD,'k-');
hold on
plot(Ns,medianSD,'r-');
legend('SD of means','SD of medians');
xlabel('Sample size');
ylabel('SD of sampling distribution');

% Formula for SD of median:
% sd_median = 1.2533 * sigma/sqrt(N)

%% Sample from a non-normal distribution
rng default;    % for reproducibility
nSamp = 15;
B = 10000;
ps = linspace(1/(B+1),1-1/(B+1),B);
tStat = zeros(B,1);
for k = 1:B
    x = datasample([-1,1],nSamp);
    tStat(k) = sqrt(nSamp) .* (mean(x) / std(x));
end
figure
tPercentiles = tinv(ps,nSamp-1);
qqplot(tStat,tPercentiles);
xlabel('Observed'); ylabel('Theoretical');