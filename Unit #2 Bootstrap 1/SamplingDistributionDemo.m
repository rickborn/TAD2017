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

figure('Name','Heights of 99 five-year-old British boys in cm');
subplot(2,2,1)
hist(H,2);
xlabel('Height (cm)');
ylabel('#');
subplot(2,2,2)
hist(H,6);
xlabel('Height (cm)');
ylabel('#');
subplot(2,2,3)
hist(H,15);
xlabel('Height (cm)');
ylabel('#');
subplot(2,2,4)
hist(H,50);
xlabel('Height (cm)');
ylabel('#');

% Jarque-Bera test for normality; 1 means we reject H0 that data is
% normally distributed
figure, qqplot(H);
[h,p] = jbtest(H);
muHgt = mean(H);
sigmaHgt = std(H);

% The true beauty and power is that, now, knowing just these two numbers, I
% can tell you the probability that a given boy is in soome height range.
% E.g. What is the probability that a boy is between 100 and 120 cm?
pHgt = normcdf(120,muHgt,sigmaHgt) - normcdf(100,muHgt,sigmaHgt);

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