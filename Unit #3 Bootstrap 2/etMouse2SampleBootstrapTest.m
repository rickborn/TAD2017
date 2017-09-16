% etMouse2SampleBootstrapTest.m
%
% RTB wrote it, working through ch. 16 of E&T, Dec. 2016

% Concepts covered:
% 1. Plotting group data with boxplot
% 2. Jittering x-values to allow plotting of raw data
% 3. Bootstrap for H0: 2 distributions are the same
% 4. Bootstrap for H0: 2 means are same (Studentized bootstrap)

%% Read in and plot the data
myAlpha = 0.05;
nBoot = 10000;

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
[pValWRST,hRS,statsRS] = ranksum(rxGrp,ctrlGrp,'Tail','right');

%% Bootstrap test for this difference
% This is identical in concept to the permutation test in that we pool the
% data and draw samples without respect to the actual group to which the
% data belong. The only difference is that we draw WITH REPLACEMENT.

H0data = [rxGrp, ctrlGrp];
nTotal = length(H0data);
muDiffBoot = zeros(nBoot,1);
for k = 1:nBoot
    xStar = H0data(unidrnd(nTotal,nTotal,1));
    %muDiffPerm(k) = trimmean(xStar(1:nRx),25) - trimmean(xStar(nRx+1:end),25);
    muDiffBoot(k) = mean(xStar(1:nRx)) - mean(xStar(nRx+1:end));
end
figure, hist(muDiffBoot,30);
hold on;
xlabel('Bootstrap differences: Rx - Ctrl'); ylabel('#');
ax = axis;
line([muDiffHat,muDiffHat],[ax(3),ax(4)],'Color','y');

%% Calculate a p-value
% This is just the # of permuted values that are greater-than-or-equal-to
% the one we observed experimentally divided by the total. E&T refer to
% this value as the "achieved significance level" or ASL.
pValBoot = sum(muDiffBoot >= muDiffHat) / nBoot;

% add this to the figure title
tStr = sprintf('E&T Fig. 16.1, p-val = %0.3f',pValBoot);
title(tStr);

%% Bootstrap to test equality of means
% The above bootstrap, as well as the permuation, examples test for the
% most radical H0, which is that both distributions are the same. But we
% might want to test the H0 that, e.g. only the means are the same. We
% cannot do this with the permutation test, but we can with the bootstrap.
% To do this, we transform each population so that they have the SAME MEAN,
% by simply subtracting each population's own mean (now zero mean), then
% adding back the pooled mean. Then we resample with replacement from each
% population and compare the t-statistic (which normalizes with the pooled
% estimate of the variance). See p. 224 of E&T

% Does it matter whether or not we add back the group mean? I get the same
% p-value in either case, and this makes sense to me.
tRxGrp = (rxGrp - mean(rxGrp)) + mean(H0data);
%tRxGrp = rxGrp - mean(rxGrp);
tCtrlGrp = (ctrlGrp - mean(ctrlGrp)) + mean(H0data);
%tCtrlGrp = ctrlGrp - mean(ctrlGrp);
tStatReal = (mean(rxGrp) - mean(ctrlGrp)) / ...
    sqrt(var(rxGrp)/nRx + var(ctrlGrp)/nCtrl);

% NOTE: E&T give the observed value of the t-statistic in the right-hand
% panel of fig. 16.1 as 0.416, but this can't be right. I get 1.0587, which
% agrees with the MATLAB t-statistic for unequal variances. i.e.:
% [hT,pValTtest,ci,statsT] = ttest2(rxGrp,ctrlGrp,'Tail','right','Vartype','unequal');
% statsT.tstat = 1.0587. In addition, my calucluated p-value is the same as
% that of E&T (~0.15).
tStatBoot = zeros(nBoot,1);
for k = 1:nBoot
    rxStar = tRxGrp(unidrnd(nRx,nRx,1));
    ctrlStar = tCtrlGrp(unidrnd(nCtrl,nCtrl,1));
    tStatBoot(k) = (mean(rxStar) - mean(ctrlStar)) / ...
       sqrt(var(rxStar)/nRx + var(ctrlStar)/nCtrl);
end
figure, hist(tStatBoot,30);
hold on;
xlabel('Studentized bootstrap differences: Rx - Ctrl'); ylabel('#');
ax = axis;
line([tStatReal,tStatReal],[ax(3),ax(4)],'Color','y');

pValTstat = sum(tStatBoot >= tStatReal) / nBoot;

% add this to the figure title
tStr = sprintf('E&T Fig. 16.1, p-val = %0.3f',pValTstat);
title(tStr);
