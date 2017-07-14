% etMouse2SampleProblem.m
%
% RTB wrote it, November 2016
% RTB adapted it for TaD 2017 on 12 July 2017

% Concepts covered:
% 1. 'boxplot' for summarizing data
% 2. 2-sample t-test with 'ttest2'
% 3. bootstrap for calculating standard error
% 4. bootstrap, parametric confidence intervals
% 5. bootstrap, percentile confidence intervals
% 6. power calculation using 'sampsizepwr'
% 7. power calculation by simulation
% 8. power curve as a function of n

%% Read in and plot the data
myAlpha = 0.05;
nBoot = 100000;

% Twenty SOD1 mice (ALS model) were randomly assigned to a treatment group
% (riluzole) or a control group. Values are their survival times, in days,
% from the beginning of treatment. 

% Read in the data
ds = dataset('xlsfile','mouseSurvivalDrugStudy');

muCtrl = mean(ds.Ctrl);
sigmaCtrl = std(ds.Ctrl);
nRx = length(ds.Rx); nCtrl = length(ds.Ctrl);

muDiffHat = mean(ds.Rx) - mean(ds.Ctrl);
% mean difference = 30.63

% plot the data
grpLabels = {'Rx';'Rx';'Rx';'Rx';'Rx';'Rx';'Rx';'Rx';'Rx';'Rx'; ...
    'Ctrl';'Ctrl';'Ctrl';'Ctrl';'Ctrl';'Ctrl';'Ctrl';'Ctrl';'Ctrl';'Ctrl'};
boxplot([ds.Rx';ds.Ctrl'],grpLabels);
xlabel('Treatment Group');
ylabel('Survival Time (days)');
tStr = sprintf('Mean Diff = %.2f days',muDiffHat);
text(1.5,180,tStr);
%boxplot([[ds.Rx,NaN,NaN]',ds.Ctrl']);
%distributionPlot([[ds.Rx,NaN,NaN]',ds.Ctrl']);

%% Read in population data
allMouseSurvival = xlsread('mouseSurvivalPopulationData.xlsx');
%% Sample from the true population:
nSamp = 10000;
muDiffSamp = zeros(nSamp,1);
for k = 1:nSamp
    rxSamp = datasample(allMouseSurvival,nRx);
    ctrlSamp = datasample(allMouseSurvival,nCtrl);
    muDiffSamp(k) = mean(rxSamp) - mean(ctrlSamp);
    
%     if muDiffSamp(k) > 30 && muDiffSamp(k) < 31
%         break;
%     end
end
figure, hist(muDiffSamp,30);
hold on
ax = axis;
hl = line([muDiffHat,muDiffHat],[ax(3),ax(4)]);
set(hl,'Color','r');
%% Use classical 2-sample t-test
[hT,pT,ci,statsT] = ttest2(ds.Rx,ds.Ctrl);
[pRS,hRS,statsRS] = ranksum(ds.Rx,ds.Ctrl);

%% Bootstrap std. error for this difference
muDiffBoot = zeros(nBoot,1);
for k = 1:nBoot
    muDiffBoot(k) = mean(ds.Rx(unidrnd(nRx,nRx,1))) - ...
        mean(ds.Ctrl(unidrnd(nCtrl,nCtrl,1)));
end
se = std(muDiffBoot);
figure, hist(muDiffBoot,30);
hold on;
xlabel('Bootstrap difference: Rx - Ctrl'); ylabel('#');
ax = axis;

%% Parametric confidence intervals:
dev = se * norminv(1-myAlpha/2);
ciParam = [muDiffHat - dev, muDiffHat + dev];
% Superimpose these on our histogram
h1 = line([ciParam(1),ciParam(1)],[ax(3),ax(4)]);
set(h1,'Color','r');
h1 = line([ciParam(2),ciParam(2)],[ax(3),ax(4)]);
set(h1,'Color','r');

%% Brute force confidence intervals:
muDiffBootSorted = sort(muDiffBoot);
iLo = ceil((myAlpha/2) * nBoot);   % index corresponding to lower bound
iHi = nBoot - iLo;                  % index corresponding to upper bound
ciBrute = [muDiffBootSorted(iLo),muDiffBootSorted(iHi)];
% Superimpose these on our histogram
h2 = line([ciBrute(1),ciBrute(1)],[ax(3),ax(4)]);
set(h2,'Color','g');
h2 = line([ciBrute(2),ciBrute(2)],[ax(3),ax(4)]);
set(h2,'Color','g');
pStr = sprintf('Parametric %d%% CI',100*(1-myAlpha));
bStr = sprintf('Brute Force %d%% CI',100*(1-myAlpha));
legend([h1,h2],pStr,bStr,'Location','northwest');

%% Added bonus demo: power calculation
% Assume mean and sigma of controls (56.22,42.48) and diff. of 30
% calculate power of this experiment at alpha = 0.05
pwrOut = sampsizepwr('t2',[56.22 42.48],86.86,[],9);    % power = 0.30

% What sample size (in each group) would we need to have power = 0.80?
nOut = sampsizepwr('t2',[56.22 42.48],86.86,0.80);      % n = 32
% 90% power?
nOut = sampsizepwr('t2',[56.22 42.48],86.86,0.90);      % n = 42

%% But note that we had to make assumptions about equal n's and sigma's. We
% can calculate power directly by simulation, using all the values
nSim = 10000;
hSim = zeros(nSim,1);
ctrlMean = mean(ds.Ctrl); ctrlSigma = std(ds.Ctrl);
rxMean = mean(ds.Rx); rxSigma = std(ds.Rx);
for k = 1:nSim
    rxSim = normrnd(rxMean,rxSigma,nRx,1);
    ctrlSim = normrnd(ctrlMean,ctrlSigma,nCtrl,1);
    [hSim(k),~] = ttest2(rxSim,ctrlSim,'Tail','right','Alpha',myAlpha);
end
myPower = sum(hSim) / nSim; % myPower = 0.2963 (with nSim = 10,000)

%% We can't directly calculate n for a given power, but we can run a
% simulation to generate a curve of power vs. n
nEach = [10:10:60];
allPower = zeros(1,length(nEach));
nSim = 1000;
for k = 1:length(nEach)
    hSim = zeros(nSim,1);
    for j = 1:nSim
        rxSim = normrnd(rxMean,rxSigma,nEach(k),1);
        ctrlSim = normrnd(ctrlMean,ctrlSigma,nEach(k),1);
        [hSim(j),~] = ttest2(rxSim,ctrlSim,'Tail','right','Alpha',myAlpha);
    end
    allPower(k) = sum(hSim) / nSim;
end
figure, plot(nEach,allPower,'k-');
hold on
plot(nEach,allPower,'r*');
xlabel('# per group'); ylabel('power');
tStr = sprintf('Power Curve (alpha = %.2f)',myAlpha);
title(tStr);