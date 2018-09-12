% etASAdemo_ng_withAnswers.m: 
% Bootstrap example for MI data from pp. 3-5 of Efron & Tibshirani
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
% 6. One-tailed vs. two-tailed tests.
% 7. CI with 'bootci' and an anonymous function
% 8. Making data tables with 'table'

% A study was done to see if low-dose aspirin would prevent heart attacks
% in healthy middle-aged men. The study design was optimal: controlled,
% randomized, double-blind. Subjects were randomly assigned to receive
% aspirin (ASA) or placebo. The summary statistics as they appeared in the
% NY Times on Jan. 27, 1987:
%
% aspirin group (n=11037): 104 heart attacks (MI)
% placebo group (n=11034): 189 heart attacks
%
% Scientific question: Does aspirin help to prevent heart attacks?

% What to do: Login to learning catalytics and join the session for the
% module entitled "ASA Bootstrap". You will answer a series of
% questions based on the guided programming below. Each section begins with
% a '%%'. Read through the comments and follow the instructions provided.
% In some cases you will be asked to answer a question, clearly indicated
% by 'QUESTION'. In other cases, you be asked to supply missing code,
% indicated by 'TODO'. The corresponding question in learning catalytics
% will be indicated in parentheses (e.g. Q1). If there is no 'Q#'
% accompanying a 'QUESTION' just type your answer into this script and
% discuss it with your team. Once you have supplied the required code, you
% can execute that section by mouse-clicking in that section (The block
% will turn yellow.) and then simultaneously hitting the 'ctrl' and 'enter'
% keys (PC) or 'command' and 'enter' keys (Mac).

%% Constants: these would normally be passed as arguments to a function
nBoot = 10000;
myAlpha = 0.05;

% useful numbers
nRx = 11037;    % number of patients in the treatment group (ASA)
nMIrx = 104;    % number of MIs in the treatment group
nCtrl = 11034;  % number of patients in the control group (placebo)
nMIctrl = 189;  % number of MIs in the control group
nTotal = nRx + nCtrl;

%% Calculate the actual ratio of rates of disease: an odds ratio

% This is defined as the ratio of 2 ratios: the numerator is the ratio of
% the number of subjects in the treatment group who had an MI divided by
% the number who did not have an MI. The denominator is the same, but for
% the control group.
% TODO: calculate the odds ratio for this study
orHat = (nMIrx / (nRx-nMIrx)) / (nMIctrl / (nCtrl-nMIctrl));

% QUESTION (Q1): What is your odds ratio to 4 decimal places?

%% Create a population from which to resample:

% The general approach in bootstrapping is to resample from our original
% sample. Thus far, we only have proportions, but we want to have the full
% information in the original sample.

% TODO: Create a column vector that is the size of each treatment group and
% that contains 1s for each person who had an MI and 0s for each person
% that did not: 1=had MI; 0=no MI
rxGrp = [ones(104,1);zeros(10933,1)];  % aspirin group for heart attacks
ctrlGrp = [ones(189,1);zeros(10845,1)];  % non-aspirin group for heart attacks

%% First method: with a 'for' loop
orStar = zeros(nBoot,1);    % holds each bootstrap calc. of the odds ratio
rng default
for k = 1:nBoot
    % TODO: Re-sample from each group WITH REPLACEMENT to create two new
    % samples: rxStar and ctrlStar. Then use these two bootstrap samples to
    % calculate an odds ratio and store it in orStar.
    rxStar = rxGrp(unidrnd(nRx,nRx,1));
    ctrlStar = ctrlGrp(unidrnd(nCtrl,nCtrl,1));
    orStar(k) = (sum(rxStar) / sum(~rxStar)) / (sum(ctrlStar) / sum(~ctrlStar));
end

%% Make a histogram of our bootstrapped ORs
figure
subplot(2,1,1);
hist(orStar,20);
hold on;
xlabel('OR^*'); ylabel('#');
title('Distribution of bootstrapped odds ratios');
% Draw a vertical line to see where our actual value lies within the
% distribution of bootstrap replications. Does it make sense?
ax = axis;
line([orHat,orHat],[ax(3),ax(4)],'Color','y','LineWidth',2);

%% Calculate the standard error and the confidence intervals

% QUESTION (Q2): What is the value of the standard error of the odds ratio?
semBoot = std(orStar);

% TODO: Use the percentile method to determine the 95% confidence interval.
orStarSorted = sort(orStar);
idxLo = floor((myAlpha/2) * nBoot);    % index corresponding to lower bound
idxHi = ceil((1-(myAlpha/2)) * nBoot); % index corresponding to upper bound
confInterval = [orStarSorted(idxLo),orStarSorted(idxHi)];

% QUESTION (Q3): What is the 95% CI based on your bootstrap distribution?

% QUESTION (Q4): What is the null value of our statistic? Ans. OR = 1

% QUESTION (Q5): How can we use the known null value of the odds ratio to
% perform a hypothesis test?

% QUESTION (Q6): Can we reject H0 at an alpha of 0.05?
% Ans. Yes. Our 95% CI does not include the null value of 1.

%% Plot CIs on histogram
line([confInterval(1),confInterval(1)],[ax(3),ax(4)],'Color','r');
line([confInterval(2),confInterval(2)],[ax(3),ax(4)],'Color','r');

%% Confidence intervals with 'bootci' and an anonymous function

% As we saw in week #2, we can also use a built-in function to calculate
% CIs with much less hassle (No 'for' loops!). In that case, we had to pass
% 'bootci' a handle to a function ('@corr'), but here we have no such
% built-in function to calculate an odds ratio. So what do we do? We make
% one up! We do this using what is known as an 'anonymous function'--a kind
% of on-the-fly function that we create in one line of code:
oddsratio = @(g1,g2) (sum(g1,'omitnan')/sum(1-g1,'omitnan'))...
    /(sum(g2,'omitnan')/sum(1-g2,'omitnan'));

% We have to make g1 and g2 be the same size:
ctrlGrp = [ctrlGrp;NaN;NaN;NaN];

% now use 'bootci':
ci = bootci(nBoot,{oddsratio,rxGrp,ctrlGrp},'alpha',myAlpha,'type','bca');

% The above is a good example of how a defect in MATLAB requires one to use
% some ingenuity in coding. I intially did this in the logical way, which
% is to use a sensible definition of the odds ratio that works perfectly
% well: oddsratio = @(g1,g2) (sum(g1)/sum(~g1))/(sum(g2)/sum(~g2));
% That is, if I now enter: oddsratio(rxGrp,ctrlGrp), I get the correct
% answer of 0.5458. However, when I try to use this anonymous function, I
% get an error: 'Nonscalar arguments to BOOTFUN must have the same number
% of rows.' That is, I can't have g1 and g2 be different lengths. So then I
% figure I'll just pad the shorter vector with NaN's:
% ctrlGrp = [ctrlGrp;NaN;NaN;NaN];
% So now they're both the same length, but then I get a different error: I
% can't use the tilde ('~') negation trick with NaN's in the vector. So now
% I have to change the function to use '1-g1' instead of '~g1' in the
% denominator of each odds calculation. So then it works, but it is rather
% unsatisfying to have to kluge things in this way!

%% Perform an explicit hypothesis test by modeling our OR under H0
% In this case, we will use a permutation test, where we resample
% WITHOUT replacement. The logic is that we are essentially randomly
% assigning each patient to the treatment or control group, then
% recalculating our odds ratio. Here, we are testing the most extreme version of
% H0, which is that the two distributions are the SAME.

% TODO: Perform resampling as though the patients all belonged to the
% same group (called H0data), shuffle this data, then arbitraily assign
% each patient to the treatment or control group and compute the odd ratio. 
% Store each bootstrapped odds ratio in orPerm

% Pool all the data:
H0data = [rxGrp;ctrlGrp];
% Place to store our results:
orPerm = zeros(nBoot,1);

rng default
for k = 1:nBoot
    shuffledData = H0data(randperm(nTotal));
    rxStar = shuffledData(1:nRx);
    ctrlStar = shuffledData(nRx+1:end);
    orPerm(k) = (sum(rxStar) / sum(~rxStar)) / (sum(ctrlStar) / sum(~ctrlStar));
end

% plot the distrubtion of permuted ORs
subplot(2,1,2);
hist(orPerm,20);
hold on
xlabel('Permuted ORs'); ylabel('#');
title('Distribution of ORs under H0');

% draw a line for our actual odds ratio
ax = axis;
line([orHat,orHat],[ax(3),ax(4)],'Color','r');
%% Calcualte a p-value
% This is a one-tailed test:
if orHat < 1
    pVal1t = sum(orPerm <= orHat) / nBoot;
else
    pVal1t = sum(orPerm >= orHat) / nBoot;
end

% The p-value can never be 0. The logic is that we could have found a
% significant value on our next iteration. Good teaching point.
if pVal1t == 0
    pVal1t = 1 / (nBoot+1);
end

% QUESTION (Q7): What is our one-tailed p-value for the odds ratio?

%% Calculate a 2-tailed p-value
if orHat < 1
    pVal2t = (sum(orPerm <= orHat) + sum(orPerm >= 1/orHat)) / nBoot;
else
    pVal2t = (sum(orPerm >= orHat) + sum(orPerm <= 1/orHat)) / nBoot;
end

if pVal2t == 0
    pVal2t = 1 / (nBoot+1);
end

% Note: This is a good teaching point. The students will initially
% calculate a 1-tailed p-value, which is the most intuitive thing to do.
% When they compare this to the results of Fisher's exact test, they may
% note that their p-value is about 1/2 of the FET p-value. Make them think
% about the other tail and how to find it. Key is that for an odds ratio,
% it is just 1/orHat.

% Philosophical interlude: The difference in implementation between
% 2-tailed and 1-tailed p-value is pretty clear. The philosophical
% difference, somewhat less. If you accidentally coded a 2-tailed test and
% get a p value of, say,0.06, and then remember "Oh! A 1-tailed test was
% actually more appropriate!" (and it really is in that instance, not for a
% "p-hacky" reason) and obtain p ~ 0.03, there's a sudden shift in
% perspective on the data. But it's the same data, and you're performing
% more or less the same analysis. Does this seem even remotely reasonable?
% This very subtle distinction would have a pretty heavy impact on a
% statistics-naïve researcher. It can be helpful to think about edge cases
% like this, where our arbitrary thresholding statistical procedure leads
% to binarization of the same data into two categories which are
% interpreted in very different ways, and how we should consider data of
% this variety. Is it helpful to construct a new categorization, e.g.
% "statistically significant (p small)," "unlikely to produce statistical
% significance (p biggish)" and "of uncertain relationship (p kinda
% small?)" or does that just move the problem?

%% Compare with the Fisher Exact Test for MI data

% load the data into an object of type 'table'
MIdata = table([nMIrx;nMIctrl],[nRx-nMIrx;nCtrl-nMIctrl],...
    'VariableNames',{'MI','NoMI'},'RowNames',{'ASA','NoASA'});

% TODO: Calculate a p-value and 95% confidence interval using Fisher's
% Exact Test
[hMI,pMI,statsMI] = fishertest(MIdata,'Tail','both','Alpha',0.05);

% h = 1
% p = 5.0328e-07
% stats = 
%              OddsRatio: 0.5458
%     ConfidenceInterval: [0.4290 0.6944]

%% Extra: Repeat calculations for a different data set:

% In the same study, the researchers also recorded the number of strokes in
% each group. Sample data for strokes: 1=had stroke; 0=no stroke
rxGrp = [ones(119,1);zeros(10918,1)];    % aspirin group for strokes
ctrlGrp = [ones(98,1);zeros(10936,1)];   % non-aspirin group for strokes

% useful numbers
nRx = length(rxGrp);      % number of patients in the treatment group (ASA)
nCtrl = length(ctrlGrp);  % number of patients in the control group (placebo)
nTotal = nRx + nCtrl;
nBoot = 10000;
myAlpha = 0.05;

% Odds ratio for stroke data:
orHat = (sum(rxGrp) / sum(~rxGrp)) / (sum(ctrlGrp) / sum(~ctrlGrp));

% TODO: Repeat the above analysis for the stroke data.
% Be sure to compare you p-values and confidence intervals that you obtain
% through bootstrapping to those obtained with Fisher's Exact Test.

strokeData = table([119;98],[11037-119;11034-98],...
    'VariableNames',{'Stroke','NoStroke'},'RowNames',{'ASA','NoASA'});
[hStroke,pStroke,statsStroke] = fishertest(strokeData,'Tail','both','Alpha',0.05);
% h = 0
% p = 0.1723
% stats = 
%              OddsRatio: 1.2163
%     ConfidenceInterval: [0.9297 1.5912]
