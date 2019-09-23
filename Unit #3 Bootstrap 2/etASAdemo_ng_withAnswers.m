% etASAdemo_ng_withAnswers.m: 
%
% Bootstrap example for Stroke data from pp. 3-5 of Efron & Tibshirani
% 
% RTB wrote it 29 October 2016 (derived from BS_ex1.m)
% RTB modified it 30 Jan 2017: combined etASAhypoth.m and etASAstats.m into
% one file that is modular for my stats class.
% RTB modified it to emphasize stroke data (13 September 2018)

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
% Scientific question #1: Does aspirin help to prevent heart attacks?
%
% In the same study, the researchers also recorded the number of strokes in
% each group:
%
% aspirin group (n=11037): 119 strokes
% placebo group (n=11034): 98 strokes
%
% Scientific question #2: Does aspirin increase the risk of having a
% stroke?
%
% We will start by addressing the 2nd question regarding strokes.

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
%
% NOTE: If you are using version R2018a of MATLAB, you won't be able to use
% the ctrl+enter feature, because it now checks the entire script for
% errors, rather than just the cell you are trying to execute. This is
% stupid, but we're stuck with it. What you can do instead is use the mouse
% to highlight the code you want to run, then hit the F9 key (PC) or you
% can also just copy the section and then paste it into the command line.

%% Constants: these would normally be passed as arguments to a function
nBoot = 10000;
myAlpha = 0.05;

% useful numbers for stroke data
nRx = 11037;        % number of patients in the treatment group (ASA)
nStrokeRx = 119;    % number of Strokes in the treatment group
nCtrl = 11034;      % number of patients in the control group (placebo)
nStrokeCtrl = 98;  % number of Strokes in the control group
figName = 'Stroke Data';

nTotal = nRx + nCtrl;

%% Calculate the actual ratio of rates of disease: an odds ratio

% This is defined as the ratio of 2 ratios: the numerator is the ratio of
% the number of subjects in the treatment group who had a stroke divided by
% the number who did not have a stroke. The denominator is the same, but
% for the control group.

% TODO: calculate the odds ratio for this study
orHat = (nStrokeRx / (nRx-nStrokeRx)) / (nStrokeCtrl / (nCtrl-nStrokeCtrl));

% QUESTION (Q1): What is your odds ratio to 4 decimal places?

%% Create a population from which to resample:

% The general approach in bootstrapping is to resample from our original
% sample. Thus far, we only have proportions, but we want to have the full
% information in the original sample.

% TODO: Create a column vector that is the size of each treatment group and
% that contains 1s for each person who had a stroke and 0s for each person
% that did not: 1=had stroke; 0=no stroke
rxGrp = [ones(nStrokeRx,1);zeros(nRx-nStrokeRx,1)];  % aspirin group for strokes
ctrlGrp = [ones(nStrokeCtrl,1);zeros(nCtrl-nStrokeCtrl,1)];  % non-aspirin group for strokes

%% Generate bootstrap replicates of the odds ratio

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

%% Make a histogram of our bootstrap replicates of OR

figure('Name',figName);
subplot(2,1,1);
hist(orStar,20);
hold on;
xlabel('OR^*'); ylabel('#');
title('Distribution of bootstrapped odds ratios');
% Draw a vertical line to see where our actual value lies within the
% distribution of bootstrap replications. Does it make sense?
ax = axis;
h1=line([orHat,orHat],[ax(3),ax(4)],'Color','g','LineWidth',2);

%% Calculate the standard error and the confidence intervals

% QUESTION (Q2): What is bootstrap estimate of the standard error of the
% odds ratio?
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
h2=line([confInterval(1),confInterval(1)],[ax(3),ax(4)],'Color','r');
line([confInterval(2),confInterval(2)],[ax(3),ax(4)],'Color','r');
legend([h1,h2],{'orHat','95% CI'},'Location','NorthEast');

%% Bonus: Confidence intervals with 'bootci' and an anonymous function

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
% NOTE: Chronux has it's own 'jackknife' function that interferes with
% MATLAB's 'jackknife' called by 'bootci' to do BCA. Solution is to remove
% it from our path for this analysis.
a = which('jackknife');
if ~contains(a,'MATLAB')    % Chronux version is interfering
    rmpath('C:\usr\rick\mat\chronux_2_11\spectral_analysis\helper');
    ci = bootci(nBoot,{oddsratio,rxGrp,ctrlGrp},'alpha',myAlpha,'type','bca');
    addpath('C:\usr\rick\mat\chronux_2_11\spectral_analysis\helper');
else
    ci = bootci(nBoot,{oddsratio,rxGrp,ctrlGrp},'alpha',myAlpha,'type','bca');
end

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

%% Bonus: 2-sided p-value from the bootstrap distribution

% We can get a more exact 2-sided p-value by gradually shrinking our % CI
% until the null value falls just outside. We have already done the hard
% work of generating the bootstrap distribution and sorting it, so now we
% just loop through and find the value of alpha that puts the null value of
% 1 outside of the lower bound.

% NOTE: This cell won't work for the MI data, since our bootstrapped
% distribution won't contain the null value.
ciLow = confInterval(1);
pAlpha = myAlpha;

while ciLow < 1
    pAlpha = pAlpha + 0.001;
    idxLo = floor((pAlpha/2) * nBoot);
    ciLow = orStarSorted(idxLo);
end
% pAlpha: 0.1620; compare with Fisher Exact value of 0.1723

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

%% Plot the distrubtion of permuted ORs

subplot(2,1,2);
hist(orPerm,20);
hold on
xlabel('Permuted ORs'); ylabel('#');
title('Distribution of ORs under H0');

% draw a line for our actual odds ratio
axis(ax);
h3=line([orHat,orHat],[ax(3),ax(4)],'Color','g','LineWidth',2);
legend(h3,'orHat','Location','NorthEast');

%% Calcualte a 1-tailed p-value

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
% statistics-na�ve researcher. It can be helpful to think about edge cases
% like this, where our arbitrary thresholding statistical procedure leads
% to binarization of the same data into two categories which are
% interpreted in very different ways, and how we should consider data of
% this variety. Is it helpful to construct a new categorization, e.g.
% "statistically significant (p small)," "unlikely to produce statistical
% significance (p biggish)" and "of uncertain relationship (p kinda
% small?)" or does that just move the problem?

%% Compare with the Fisher Exact Test for Stroke data

strokeData = table([nStrokeRx;nStrokeCtrl],[nRx-nStrokeRx;nCtrl-nStrokeCtrl],...
    'VariableNames',{'Stroke','NoStroke'},'RowNames',{'ASA','NoASA'});
[hStroke,pStroke,statsStroke] = fishertest(strokeData,'Tail','both','Alpha',0.05);
% h = 0
% p = 0.1723
% stats = 
%              OddsRatio: 1.2163
%     ConfidenceInterval: [0.9297 1.5912]

%% Compare with a Chi-square test

observedData = [nStrokeRx,nRx-nStrokeRx;nStrokeCtrl,nCtrl-nStrokeCtrl];

% Calculate our expected table based on row and column totals:
CT = sum(observedData,1);   % column totals
RT = sum(observedData,2);   % row totals
GT = sum(RT);               % grand total
expectedData = (repmat(CT,length(RT),1) ./ GT) .* repmat(RT,1,length(CT));

% Calculate our chi2 statistic:
chi2stat = sum(((expectedData(:) - observedData(:)).^2) ./ expectedData(:));

% Calculate our degrees of freedom:
df = (size(observedData,1)-1) * (size(observedData,2)-1);
% or, more cleverly:
%df = prod(size(observedData)-1);

% Use MATLAB's built-in chi2 distribution to get a p-value:
pChi2 = chi2cdf(chi2stat,df,'upper');

%% Repeat calculations for the heart attack data:

% useful numbers for heart attack data
nRx = 11037;    % number of patients in the treatment group (ASA)
nMIrx = 104;    % number of MIs in the treatment group
nCtrl = 11034;  % number of patients in the control group (placebo)
nMIctrl = 189;  % number of MIs in the control group
nTotal = nRx + nCtrl;
figName = 'MI Data';

nBoot = 10000;
myAlpha = 0.05;

% Generate raw data
rxGrp = [ones(nMIrx,1);zeros(nRx-nMIrx,1)];         % aspirin group for MI
ctrlGrp = [ones(nMIctrl,1);zeros(nCtrl-nMIctrl,1)]; % non-aspirin group for MI

% Odds ratio for MI data:
orHat = (sum(rxGrp) / sum(~rxGrp)) / (sum(ctrlGrp) / sum(~ctrlGrp));

% TODO: Repeat the above analysis for the MI data. Go back to line 100 and
% re-run the code with the new numbers as above.
%
% Be sure to compare you p-values and confidence intervals that you obtain
% through bootstrapping to those obtained with Fisher's Exact Test.
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

%% Final thoughts
% 
% There is a beautiful and important section on the relationship between
% confidence intervals and hypothesis tests in section 15.4 of E&T (p.
% 214). 
% 
% Look at the top half of your figure. Key question: What value of alpha
% will make the lower end of the bootstrap confidence interval equal to 1?
% It is just the probability mass of the bootstrap distribution to the left
% of 1: sum(orStar <= 1) / nBoot. This gives us a bootstrap p-value of
% 0.0807. (Compare this with our 1-tailed permutation-based p-value of
% 0.0887.)
% 
% Compare the top and bottom halves and consider this statement from E&T
% (p. 216): "In this sense, the permutation p-value measures how far the
% observed estimate, orHat, is from 1, while the bootstrap p-value
% measures how far 1 is from orHat." (Remember that 1 is our null value
% for the odds ratio.)
% 
% Which test is preferred? In general, the permutation test is the most
% rigorous way to test H0, since we are directly instantiating it by
% "breaking" the relationship we seek to test. In this case, we are
% randomly assigning the observed data to either the treatment or the
% control group. As E&T point out (p. 216): "The permutation p-value is
% exact, while the bootstrap p-value is approximate." They conclude with
% what should become one of your statistical mantras (p. 218): "When there
% is something to permute, it is a good idea to do so . . . ."