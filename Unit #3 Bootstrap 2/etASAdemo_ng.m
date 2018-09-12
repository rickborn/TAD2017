% etASAdemo_ng.m: Bootstrap example for MI data from pp. 3-5 of Efron &
% Tibshirani
% 
% RTB wrote it 29 October 2016 (derived from BS_ex1.m)
% RTB modified it 30 Jan 2017: combined etASAhypoth.m and etASAstats.m 
% into one file that is modular for my stats class.

% Concepts covered: 
% 1. Test for proportions: odds ratio 
% 2. Comparing resampling tests with Fisher's exact test 
% 3. Std. error and confidence intervals through bootstrapping
% 4. Relationship between CI and hypothesis test 
% 5. Permutation test for strong test of H0. 
% 6. One-tailed vs. two-tailed tests.

% A study was done to see if low-dose aspirin would prevent heart attacks
% in healthy middle-aged men. The study design was optimal: controlled,
% randomized, double-blind. Subjects were randomly assigned to receive
% aspirin (ASA) or placebo. The summary statistics as they appeared in the
% NY Times on Jan. 27, 1987:
%
% aspirin group (n=11037): 104 heart attacks (MI) 
% placebo group (n=11034): 189 heart attacks
%
% Does aspirin help to prevent heart attacks?

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
myAlpha = 0.05;     % for a 95% confidence interval

% useful numbers
nRx = 11037;    % number of patients in the treatment group (ASA)
nMIrx = 104;    % number of MIs in the treatment group
nCtrl = 11034;  % number of patients in the control group (placebo)
nMIctrl = 189;  % number of MIs in the control group
nTotal = nRx + nCtrl;

%% Calculate the actual ratio of rates of disease: an odds ratio

% This is defined as the ratio of 2 ratios. The 'odds' of an event with a
% binary outcome (i.e. it either happens or it doesn't) is defined as the
% probability that the event happens divided by the probability that it
% does not happen. So we calculate the odds for each condition (treatment
% vs. control) and define the odds ratio as the odds of having a heart
% attack given that you were in the treatment group divided by the odds of
% having a heart attack given that you were in the control group.

% TODO: calculate the odds ratio for this study
orHat = ;

% QUESTION (Q1): What is your odds ratio to 4 decimal places?

%% Create a population from which to resample:

% The general approach in bootstrapping is to resample from our original
% sample. Thus far, we only have proportions, but we want to have the full
% information in the original sample.

% TODO: Create a column vector that is the size of each group and
% that contains 1s for each person who had an MI and 0s for each person
% who did not: 1=had MI; 0=no MI. HINT: length(rxGrp) should = 11037.
rxGrp = ;  % aspirin group for heart attacks
ctrlGrp = ;  % non-aspirin group for heart attacks

%% First method: with a 'for' loop
orStar = zeros(nBoot,1);    % holds each bootstrap calc. of the odds ratio
rng default
for k = 1:nBoot
    % TODO: Re-sample from each group WITH REPLACEMENT to create two new
    % samples: rxStar and ctrlStar. Then use these two bootstrap samples to
    % calculate an odds ratio and store it in orStar.
    rxStar = ;
    ctrlStar = ;
    orStar(k) = ;
end

%% Make a histogram of our bootstrapped ORs
figure, hist(orStar,20);
hold on;
xlabel('OR^*'); ylabel('#');
title('Distribution of bootstrapped odds ratios');
% Draw a vertical line to see where our actual value lies within the
% distribution of bootstrap replications. Does it make sense?
ax = axis;
line([orHat,orHat],[ax(3),ax(4)],'Color','y');
%% Calculate the standard error and the confidence intervals

% QUESTION (Q2): What is the value of the standard error of the odds ratio?
semBoot = ;

% TODO: Use the percentile method to determine the 95% confidence interval
% based on your distribution of bootstrap replications.

confInterval = ;

% QUESTION (Q3): What is the 95% CI based on your bootstrap distribution?

% QUESTION (Q4): What is the null value of our statistic? Ans. OR = 1

% QUESTION (Q5): How can we use the known null value of the odds ratio to
% perform a hypothesis test?

% QUESTION (Q6): Can we reject H0 at an alpha of 0.05?

%% Plot CIs on histogram
line([confInterval(1),confInterval(1)],[ax(3),ax(4)],'Color','r');
line([confInterval(2),confInterval(2)],[ax(3),ax(4)],'Color','r');

%% Perform an explicit hypothesis test by modeling our OR under H0

% In this case, we will use a permutation test, where we resample WITHOUT
% replacement. The logic is that we are essentially randomly assigning each
% patient to the treatment or control group, then recalculating our odds
% ratio. Here, we are testing the most extreme version of H0, which is that
% the two distributions are the SAME.

% TODO: Perform resampling as though the patients all belonged to the
% same group (called H0data), shuffle this data, then arbitraily assign
% each patient to the treatment or control group and compute the odds ratio. 
% Store each bootstrapped odds ratio in orPerm

% Pool all the data:
H0data = [rxGrp;ctrlGrp];
% Place to store our results:
orPerm = zeros(nBoot,1);

rng default
for k = 1:nBoot
    % Shuffle and assign to rxStar and ctrlStar (ALL values are used!)
    
    % Compute the odds ratio for the shuffled data
    orPerm(k) = ;
end

% plot the distrubtion of permuted ORs
figure, hist(orPerm,20);
hold on
xlabel('Permuted ORs'); ylabel('#');
title('Distribution of ORs under H0');

% Draw a line for our actual odds ratio
ax = axis;
line([orHat,orHat],[ax(3),ax(4)],'Color','r');
%% Calcualte a one-tailed p-value

% TODO: Calculate a one-tailed p-value based on your permuted samples (i.e.
% orPerm) and store it in a variable called 'pVal1t'


% QUESTION (Q7): What is our one-tailed p-value for the odds ratio?

%% Calculate a 2-tailed p-value

% TODO: Calculate a two-tailed p-value based on your permuted samples (i.e.
% orPerm) and store it in a variable called 'pVal2t'

% QUESTION (Q8): What is our two-tailed p-value for the odds ratio?

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

% TODO: Calculate a 2-tailed p-value and 95% confidence interval using
% Fisher's Exact Test

% QUESTION (Q9): What p-value does Fisher's Exact Test give?

% QUESTION (Q10): Why might your results differ from Fisher's Exact Test

% QUESTION (Q11): What is the 95% CI from Fisher's Exact Test?
% Be sure to compare this with your values from the bootstrap!

%% Extra: Repeat calculations for a different data set:

% In the same study, the researchers also recorded the number of strokes in
% each group. Here are the data:
%
% Treatment group:
% 119 had a stroke; 10918 did not have a stroke
%
% Control group:
% 98 had a stroke; 10936 did not have a stroke

% TODO: Repeat the above analysis for the stroke data.
% Be sure to compare the p-values and confidence intervals that you obtain
% through bootstrapping to those obtained with Fisher's Exact Test.
