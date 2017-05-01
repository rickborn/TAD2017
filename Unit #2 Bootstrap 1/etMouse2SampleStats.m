function [se,ciParam,ciBrute] = etMouse2SampleStats(nBoot,myAlpha)

% etMouse2SampleStats.m: computes standard error and conficendce intervals
%
% RTB wrote it, November 2016.
% Note: This is the function version of 'etMouse2SampleDemo.m'

% Concepts covered:
% 1. 'boxplot' for summarizing data
% 2. 2-sample t-test with 'ttest2'
% 3. bootstrap for calculating standard error
% 4. bootstrap, parametric confidence intervals
% 5. bootstrap, percentile confidence intervals
% 6. power calculation using 'sampsizepwr'
% 7. power calculation by simulation
% 8. power curve as a function of n

if nargin < 2, myAlpha = 0.05; end
if nargin < 1, nBoot = 10000; end

% Sixteen SOD1 mice (ALS model) were randomly assigned to a treatment group
% (riluzole) or a control group. Values are their survival times, in days,
% from the beginning of treatment. [RTB made this up.]
rxGrp = [94,197,16,38,99,141,23];   % mean of 86.857, sigma 66.767
ctrlGrp = [52,104,146,10,51,30,40,27,46]; % mean of 56.222, sigma 42.476

nRx = length(rxGrp); nCtrl = length(ctrlGrp);

muDiffHat = mean(rxGrp) - mean(ctrlGrp);
% mean difference = 30.63

% Bootstrap std. error for this difference
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

% Parametric confidence intervals:
dev = se * norminv(1-myAlpha/2);
ciParam = [muDiffHat - dev, muDiffHat + dev];
% Superimpose these on our histogram
h1 = line([ciParam(1),ciParam(1)],[ax(3),ax(4)]);
set(h1,'Color','r');
h1 = line([ciParam(2),ciParam(2)],[ax(3),ax(4)]);
set(h1,'Color','r');

% Brute force confidence intervals:
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

% Added bonus demo: power calculation
% Assume mean and sigma of controls (56.22,42.48) and diff. of 30
% calculate power of this experiment at alpha = 0.05
pwrOut = sampsizepwr('t2',[56.22 42.48],86.86,[],9);    % power = 0.30

% What sample size (in each group) would we need to have power = 0.80?
nOut = sampsizepwr('t2',[56.22 42.48],86.86,0.80);      % n = 32
% 90% power?
nOut = sampsizepwr('t2',[56.22 42.48],86.86,0.90);      % n = 42

% But note that we had to make assumptions about equal n's and sigma's. We
% can calculate power directly by simulation, using all the values
nSim = 10000;
hSim = zeros(nSim,1);
ctrlMean = mean(ctrlGrp); ctrlSigma = std(ctrlGrp);
rxMean = mean(rxGrp); rxSigma = std(rxGrp);
for k = 1:nSim
    rxSim = normrnd(rxMean,rxSigma,nRx,1);
    ctrlSim = normrnd(ctrlMean,ctrlSigma,nCtrl,1);
    [hSim(k),~] = ttest2(rxSim,ctrlSim,'Tail','right','Alpha',myAlpha);
end
myPower = sum(hSim) / nSim; % myPower = 0.2963 (with nSim = 10,000)

% We can't directly calculate n for a given power, but we can run a
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