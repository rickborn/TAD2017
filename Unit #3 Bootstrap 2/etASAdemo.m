% etASAdemo.m: Bootstrap example for stroke data from pp. 3-5 of Efron & Tibshirani
% 
% RTB wrote it 29 October 2016 (derived from BS_ex1.m)
% RTB modified it 30 Jan 2017: combined etASAhypoth.m and etASAstats.m into
% one file that is modular for my stats class.

% Concepts covered:
% 1. Test for proportions: odds ratio
% 2. Comparing resampling tests with Fisher's exact test
% 3. Std. error and confidence intervals through bootstrapping
% 4. Relationship between CI and hypothesis test
% 5. Permutation test for strong test of H0.
% 6. One-sided vs. two-sided tests.

% The basic idea is simple. For each study we have a large number of
% patients in each group: Those treated with aspririn (=ASA) and those
% treated with placebo. These patients were folowed for some length of
% time, and the number within each group that either had or did not have a
% stroke (or MI) were counted. For the bootstrap, we create an array for
% each group of patients, with '1' meaning the patient had a stroke/MI and
% '0' meaning the patient did not. For the boot strap test, we then
% randomly sample patients WITH REPLACEMENT to get a new, so-called 
% 'bootstrap sample'. We do this a bunch of times, then get a distribution
% of descriptive statistics on the boot strap sample. In this case, we are
% simply asking about the ratio of the disease rates in the two groups.

% Constants: these would normally be passed as arguments to a function
nBoot = 10000;
myAlpha = 0.05;

%% We will just use the MI data for this example.
% sample data for heart attacks (MI): 1=had MI; 0=no MI
rxGrp = [ones(104,1);zeros(10933,1)];  % aspirin group for heart attacks
ctrlGrp = [ones(189,1);zeros(10845,1)];  % non-aspirin group for heart attacks

% sample data for strokes: 1=had stroke; 0=no stroke
% rxGrp = [ones(119,1);zeros(10918,1)];  % aspirin group for strokes
% ctrlGrp = [ones(98,1);zeros(10936,1)];   % non-aspirin group for strokes

nRx = length(rxGrp);
nCtrl = length(ctrlGrp);
nTotal = nRx + nCtrl;

%% Calculate the actual ratio of rates of disease: an odds ratio
% Note that in the original E&T version, the statistic is computed as the
% number having the event divided by the total. But the more traditional
% way (for comparison with Fisher's exact test) is to put the number not
% having the event in the denominator.
orHat = (sum(rxGrp) / sum(~rxGrp)) / (sum(ctrlGrp) / sum(~ctrlGrp));

%% First method: with a 'for' loop
orStar = zeros(nBoot,1);    % holds each bootstrap calc. of the odds ratio

for k = 1:nBoot
    % Grab two samples of the appropriate size from the pooled data
    rxStar = rxGrp(unidrnd(nRx,nRx,1));
    ctrlStar = ctrlGrp(unidrnd(nCtrl,nCtrl,1));
    orStar(k) = (sum(rxStar) / sum(~rxStar)) / (sum(ctrlStar) / sum(~ctrlStar));
end

%% Make a histogram of our bootstrapped ORs
figure
hist(orStar,20);
hold on;
xlabel('OR^*'); ylabel('#');
title('Distribution of bootstrapped odds ratios');
% draw a vertical line for our actual value
ax = axis;
line([orHat,orHat],[ax(3),ax(4)],'Color','y');
%% Calculate the standard error and the confidence intervals
stdError = std(orStar);

% determine confidence interval, CI
C = sort(orStar);
idxLo = floor((myAlpha/2) * nBoot);   % index corresponding to lower bound
idxHi = nBoot - idxLo;                  % index corresponding to upper bound
confInterval = [C(idxLo),C(idxHi)];
% Q: What is the null value of our statistic? (Ans: OR = 1)
% Q: Can we use this to perform a hypothesis test? (Ans: Yes. Since our 95%
% CI does not include the null value of 1, we can reject H0 at alpha=0.05.
%% Plot CIs on histogram
line([confInterval(1),confInterval(1)],[ax(3),ax(4)],'Color','r');
line([confInterval(2),confInterval(2)],[ax(3),ax(4)],'Color','r');

%% Perform an explicit hypothesis test by modeling our OR under H0
% In this case, we prefer to use a permutation test, where we resample
% WITHOUT replacement. The logic is that we are essentially randomly
% assigning each patient to the treatment or control group, then
% recalculating our OR. Here, we are testing the most extreme version of
% H0, which is that the two distributions are the same.
H0data = [rxGrp;ctrlGrp];
orPerm = zeros(nBoot,1);
for k = 1:nBoot
    shuffledData = H0data(randperm(nTotal));
    rxStar = shuffledData(1:nRx);
    ctrlStar = shuffledData(nRx+1:end);
    orPerm(k) = (sum(rxStar) / sum(~rxStar)) / (sum(ctrlStar) / sum(~ctrlStar));
end

% plot the distrubtion of permuted ORs
figure
hist(orPerm,20);
hold on
xlabel('Permuted ORs'); ylabel('#');
title('Distribution of ORs under H0');
ax = axis;
line([orHat,orHat],[ax(3),ax(4)],'Color','r');
%% Calcualte a p-value
% This is a one-sided test:
if orHat < 1
    pVal1s = sum(orPerm <= orHat) / nBoot;
else
    pVal1s = sum(orPerm >= orHat) / nBoot;
end

% The p-value can never be 0. The logic is that we could have found a
% significant value on our next iteration. Good teaching point.
if pVal1s == 0
    pVal1s = 1 / (nBoot+1);
end

%% Calculate a 2-sided p-value
if orHat < 1
    pVal2s = (sum(orPerm <= orHat) + sum(orPerm >= 1/orHat)) / nBoot;
else
    pVal2s = (sum(orPerm >= orHat) + sum(orPerm <= 1/orHat)) / nBoot;
end

if pVal2s == 0
    pVal2s = 1 / (nBoot+1);
end

% Note: This is a good teaching point. The students will initially
% calculate a 1-sided p-value, which is the most intuitive thing to do.
% When they compare this to the results of Fisher's exact test, they may
% note that their p-value is about 1/2 of the FET p-value. Make them think
% about the other tail and how to find it. Key is that for an odds ratio,
% it is just 1/orHat.
%% Compare with the Fisher Exact Test for MI data
MIdata = table([104;189],[11037-104;11034-189],...
    'VariableNames',{'MI','NoMI'},'RowNames',{'ASA','NoASA'});
[hMI,pMI,statsMI] = fishertest(MIdata,'Tail','both','Alpha',0.05);
% h =
%      1
% p =
%    5.0328e-07
% stats = 
%              OddsRatio: 0.5458
%     ConfidenceInterval: [0.4290 0.6944]

%% Compare with the Fisher Exact Test for stroke data
strokeData = table([119;98],[11037-119;11034-98],...
    'VariableNames',{'Stroke','NoStroke'},'RowNames',{'ASA','NoASA'});
[hStroke,pStroke,statsStroke] = fishertest(strokeData,'Tail','both','Alpha',0.05);
% h =
%      0
% p =
%     0.1723
% stats = 
%              OddsRatio: 1.2163
%     ConfidenceInterval: [0.9297 1.5912]