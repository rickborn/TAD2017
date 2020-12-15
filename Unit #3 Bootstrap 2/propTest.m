function [tHat,tCI,pVal,dT] = propTest(nRx,nErx,nCtrl,nEctrl,statFlag,nBoot,myAlpha,pFlag)

% propTest.m: bootstrap test of proportions, compare with Fisher Exact Test
%
% [tHat,tCI,pVal,dT] = propTest(nRx,nErx,nCtrl,nEctrl,statFlag,nBoot,myAlpha,pFlag)
% e.g. [tHat,tCI,pVal,dT] = propTest(1000,8,1000,86,1,10000,0.05,1);
%
% Inputs:
% - nRx, # of subjects in the treatment group
% - nErx, # of events in the treatment group
% - nCrtl, # of subjects in the control group
% - nEctrl, # of events in the control group
% - statFlag, 1 = 'odds ratio' (default), 0 = 'efficacy' as test statistic
% - nBoot, # of bootstrap iterations (default = 10,000)
% - myAlpha, alpha value for confidence interval (default = 0.05)
% - pFlag, 1 = plot histograms (default = 0)
%
% Outputs:
% - tHat, value of test statistic computed from data
% - tCI, confidence interval for test statistic (default = 95% CI)
% - pVal, p-value for permutation test
% - dT, data table, formatted for Fisher Exact Test (fishertest.m)
%
% RTB wrote it, 09 November 2020; itchy hives after Nahant swim
% added option to use 'efficacy' as test statistic, because this is what
% the vaccine people use. Numbers from Pfizer and Moderna look good!
%
% modified from etASAdemo_ng_withAnswers.m

%% Defaults
if nargin < 8, pFlag = 0; end
if nargin < 7, myAlpha = 0.05; end
if nargin < 6, nBoot = 10000; end
if nargin < 5, statFlag = 1; end
if nargin < 4, error('Data needed!'); end
    
%% Set up numbers for bootstrap

nTotal = nRx + nCtrl;
rxGrp = [ones(nErx,1);zeros(nRx-nErx,1)];           % Events in treatment group
ctrlGrp = [ones(nEctrl,1);zeros(nCtrl-nEctrl,1)];   % Events in control group

% Format into a table for the Fisher Exact Test
dT = table([nErx;nEctrl],[nRx-nErx;nCtrl-nEctrl],...
    'VariableNames',{'Event','No Event'},'RowNames',{'Rx','Ctrl'});

%% Determine which test statistic to use:

% Ultimately we may want to use a switch/case statement to allow for a
% variety of different test statistics. Right now we have 2 options: the
% 'odds ratio' or the 'efficacy'
if statFlag
    tStat = @(g1,g2) (sum(g1)/sum(~g1))/(sum(g2)/sum(~g2)); % odds ratio
    tName = 'odds ratio';
else
    tStat = @(g1,g2) 1 - ((sum(g1)/length(g1)) / (sum(g2)/length(g2))); % efficacy
    tName = 'efficacy';
end

%% Calculate the actual statistic from our data
tHat = tStat(rxGrp,ctrlGrp);

%% Generate bootstrap replicates of the test statistic

tStar = zeros(nBoot,1); % holds each bootstrap calc. of the test statistic
for k = 1:nBoot
    rxStar = rxGrp(unidrnd(nRx,nRx,1));
    ctrlStar = ctrlGrp(unidrnd(nCtrl,nCtrl,1));
    tStar(k) = tStat(rxStar,ctrlStar);
end

%% Make a histogram of our bootstrap replicates of OR

if pFlag
    figure;
    subplot(2,1,1);
    histogram(tStar,20);
    hold on;
    xlabel('t^*'); ylabel('#');
    title(['Distribution of bootstrapped ',tName]);
    % Draw a vertical line to see where our actual value lies within the
    % distribution of bootstrap replications. Does it make sense?
    ax = axis;
    h1=line([tHat,tHat],[ax(3),ax(4)],'Color','g','LineWidth',2);
end

%% Use the percentile method to determine the 95% confidence interval.

tStarSorted = sort(tStar);
idxLo = floor((myAlpha/2) * nBoot);    % index corresponding to lower bound
idxHi = ceil((1-(myAlpha/2)) * nBoot); % index corresponding to upper bound
tCI = [tStarSorted(idxLo),tStarSorted(idxHi)];

%% Plot CIs on histogram
if pFlag
    h2=line([tCI(1),tCI(1)],[ax(3),ax(4)],'Color','r');
    line([tCI(2),tCI(2)],[ax(3),ax(4)],'Color','r');
    legend([h1,h2],{'tHat','95% CI'},'Location','NorthEast');
end

%% Perform an explicit hypothesis test by modeling our OR under H0

% In this case, we will use a permutation test, where we resample WITHOUT
% replacement. The logic is that we are essentially randomly assigning each
% patient to the treatment or control group, then recalculating our odds
% ratio. Here, we are testing the most extreme version of H0, which is that
% the two distributions are the SAME.

% Pool all the data:
H0data = [rxGrp;ctrlGrp];
% Place to store our results:
tPerm = zeros(nBoot,1);

rng default
for k = 1:nBoot
    shuffledData = H0data(randperm(nTotal));
    rxStar = shuffledData(1:nRx);
    ctrlStar = shuffledData(nRx+1:end);
    tPerm(k) = tStat(rxStar,ctrlStar);
end

%% Plot the distrubtion of permuted ORs

if pFlag
    subplot(2,1,2);
    histogram(tPerm,20);
    hold on
    xlabel(['Permuted ' tName]); ylabel('#');
    title(['Distribution of ' tName ' under H0']);
    
    ax = axis;
    h3=line([tHat,tHat],[ax(3),ax(4)],'Color','g','LineWidth',2);
    legend(h3,'tHat','Location','NorthEast');
end

%% Calculate a 2-tailed p-value

% The null value for the odds ratio is 1
if statFlag
    if tHat < 1
        pVal2t = (sum(tPerm <= tHat) + sum(tPerm >= 1/tHat)) / nBoot;
    else
        pVal2t = (sum(tPerm >= tHat) + sum(tPerm <= 1/tHat)) / nBoot;
    end
% The null value for the efficacy is 0
else
    if tHat > 0
        pVal2t = (sum(tPerm >= tHat) + sum(tPerm <= -tHat)) / nBoot;
    else
        pVal2t = (sum(tPerm <= tHat) + sum(tPerm >= -tHat)) / nBoot;
    end
end

% p-value can never be 0
if pVal2t == 0
    pVal2t = 1 / (nBoot+1);
end
pVal = pVal2t;