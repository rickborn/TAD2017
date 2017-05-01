% etMouse2SamplePermutationTest.m
%
% RTB wrote it, working through Ch. 15 of E&T, Dec. 2017
%
% E&T Wisdom: "When there is something to permute, it is a good idea to do
% so . . . ."

% Concepts covered:
% 1. Using 'boxplot' to visualize grouped data
% 2. Jittering x-values to allow visualization of raw data
% 3. Permutation test for hypothesis testing
% 4. Comparison with 2-sample t-test
% 5. Bootstrap test for hypothesis testing
% 6. Permutation vs. bootstrap

%% Read in and plot the data
myAlpha = 0.05;
nPerm = 10000;

% Sixteen SOD1 mice (ALS model) were randomly assigned to a treatment group
% (riluzole) or a control group. Values are their survival times, in days,
% from the beginning of treatment. [RTB made this up.]
rxGrp = [94,197,16,38,99,141,23];   % mean of 86.857, sigma 66.767
ctrlGrp = [52,104,146,10,51,30,40,27,46]; % mean of 56.222, sigma 42.476

nRx = length(rxGrp); nCtrl = length(ctrlGrp);

muDiffHat = mean(rxGrp) - mean(ctrlGrp);
% mean difference = 30.63

% plot the data
grpLabels = {'Rx';'Rx';'Rx';'Rx';'Rx';'Rx';'Rx'; ...
    'Ctrl';'Ctrl';'Ctrl';'Ctrl';'Ctrl';'Ctrl';'Ctrl';'Ctrl';'Ctrl'};
boxplot([rxGrp';ctrlGrp'],grpLabels);
hold on
xlabel('Treatment Group');
ylabel('Survival Time (days)');
tStr = sprintf('Mean Diff = %.2f days',muDiffHat);
text(1.5,180,tStr);
%boxplot([[rxGrp,NaN,NaN]',ctrlGrp']);
%figure, distributionPlot([[rxGrp,NaN,NaN]',ctrlGrp']);

%% Superimpose raw data
% To do this, we need to know the numeric values associated with the two
% groups. This can be obtained using the 'gca' command
xVals = get(gca,'XTick');

% It might help to jitter things a bit:
jFactor = (xVals(2) - xVals(1)) / 50;
plot((ones(1,nRx) .* xVals(1)) + (randn(1,nRx).*jFactor),rxGrp,'ko');
plot((ones(1,nCtrl) .* xVals(2)) + (randn(1,nCtrl).*jFactor),ctrlGrp,'ko');

%% Use classical 2-sample t-test and a non-parametric test
[hT,pValTtest,ci,statsT] = ttest2(rxGrp,ctrlGrp,'Tail','right');
% [pValWRST,hRS,statsRS] = ranksum(rxGrp,ctrlGrp,'Tail','right');

%% Permutation test for this difference
% Key step: Our null hypothesis is that the two data sets came from the
% SAME underlying probability distribution. Put another way, "if H0 is true,
% then any of our observations (survival times) for any of the mice could
% have come equally well from either of the treatment groups." (E&T p.
% 205). So we just pool the data, shuffle it, then take the first nRx
% values to be our Rx group (under H0) and the remaining to be the ctrl
% group.Note that this works for any statistic. For example, compare
% looking at the difference with 'mean' vs. 'trimmean' (the trimmed mean).
H0data = [rxGrp, ctrlGrp];
nTotal = length(H0data);
muDiffPerm = zeros(nPerm,1);
for k = 1:nPerm
    shuffledData = H0data(randperm(nTotal));
    %muDiffPerm(k) = trimmean(shuffledData(1:nRx),25) - trimmean(shuffledData(nRx+1:end),25);
    muDiffPerm(k) = mean(shuffledData(1:nRx)) - mean(shuffledData(nRx+1:end));
end
figure, hist(muDiffPerm,30);
hold on;
xlabel('Permuted differences: Rx - Ctrl'); ylabel('#');
ax = axis;
line([muDiffHat,muDiffHat],[ax(3),ax(4)],'Color','y');

%% Calculate a p-value
% This is just the # of permuted values that are greater-than-or-equal-to
% the one we observed experimentally divided by the total. E&T refer to
% this value as the "achieved significance level" or ASL.
pValPerm = sum(muDiffPerm >= muDiffHat) / nPerm;

% add this to the figure title
tStr = sprintf('E&T Fig. 15.1, p-val = %0.3f',pValPerm);
title(tStr);

%% Use bootstrap to do the hypothesis test
% Here the strategy is to resample from each group separately, then look at
% how the difference behaves.
nBoot = 10000;
muDiffBoot = zeros(nBoot,1);
for k = 1:nBoot
    muDiffBoot(k) = mean(rxGrp(unidrnd(nRx,nRx,1))) - ...
        mean(ctrlGrp(unidrnd(nCtrl,nCtrl,1)));
end
se = std(muDiffBoot);
figure, hist(muDiffBoot,30);
hold on;
xlabel('Bootstrap difference: Rx - Ctrl'); ylabel('#');
ax = axis;
line([muDiffHat,muDiffHat],[ax(3),ax(4)],'Color','y');
line([0,0],[ax(3),ax(4)],'Color','r');

pValBS = sum(muDiffBoot < 0) / nBoot;
tStr = sprintf('E&T Fig. 15.4, p-val = %0.3f',pValBS);
title(tStr);

% There is a beautiful and important section on the relationship between
% confidence intervals and hypothesis tests in section 15.4 of E&T (p.
% 214). Key question: What value of alpha will make the lower end of the
% bootstrap confidence interval equal to 0? It is just the probability mass
% of the bootstrap distribution to the left of 0:
% sum(muDiffBoot < 0) / nBoot
%
% And another gem on p. 216: "In this sense, the permutation p-value
% measures how far the observed estimate, theta_hat, is from 0, while the
% bootstrap p-value measures how far 0 is from theta_hat.
%
% And another (p. 216): "The permutation p-value is exact, while the
% bootstrap p-value is approximate."
%
% And, finally, the gem of gems (p. 218): "When there is something to permute, it is
% a good idea to do so . . . ."
