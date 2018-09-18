% FWERwithCorrSim.m
%
% solution code for Thornquist test questions
% RTB wrote it, 27 Oct 2017

% Questions 20-28 will be related to the same data ("anxiety.mat" in the
% Desire2Learn Week 8 folder). Each patient is represented as a row, and
% their normalized score on each test occupies each column (i.e. row 1 is
% Patient 1, and the first column of row 1 is Patient 1's score on Test 1).
% 
% You are a researcher studying anxiety in humans, and you want to know if
% a particular mutation makes people more likely to have anxiety. Each
% patient takes 5 different tests that measure anxiety: 20 patients with a
% known mutation ("experimental"), and 20 patients with the wild type
% allele ("control"). You then normalize the data so that a score of 0 is
% the population mean on the test, and a score of 1 is one standard
% deviation away from the population mean, and use two-sample t-tests to
% discern if there is a difference between the means of the two groups on
% each test. The experimental and control scores are plotted here for each
% test (mean +/- S.E.M.).

%% Load in data

load anxiety_assays.mat

%% Plot the data

% anonymous function for standard error of the mean
sem = @(x) std(x) ./ sqrt(length(x));

figure('Name','Anxiety Test Scores');
h1=bar(mean(expt),'b');
hold on
errorbar(mean(expt),sem(expt),'b.');
h2=bar(mean(ctl),'r');
errorbar(mean(ctl),sem(ctl),'r.');
legend([h1,h2],{'Experimental','Control'});
xlabel('Test #'); ylabel('Score');

%% Do separate t-tests for each test

% Unadjusted p-values:
[H,pVals] = ttest2(expt,ctl);

% Q15: Without making any post-hoc adjustments for multiple comparisons, on
% which tests is there a significant difference between the "experimental"
% and "control" treatments (at a significance level of p < 0.05)? (Check
% all that apply)
%
% A15: All p-values are below 0.05

% Q16:  Using the standard Bonferroni correction are all of the tests still
% significant at p < 0.05? Ans. No
% A16: We can adjust things in one of two ways:
%   1. divide our criterion by the number of tests and use the same
%      p-values
%   2. multiply our p-values by the number of tests and use the same
%      criterion
nTests = length(pVals);
BCpVals = pVals .* nTests;

% A16: No. Test #3 now gives a corrected p-value of 0.0578.

% Q17: Should you be concerned about multiple comparisons here?
% A17: Yes; we performed multiple hypothesis tests, and so we have a higher
% probability of at least one test rejecting the null hypothesis, even if
% it's true, and so we should use a more stringent criterion.

%% Assess correlation among tests:

% The purpose of multiple comparisons corrections is to address the fact
% that, if you perform several independent hypothesis tests, each at a
% particular level of significance, the probability that at least one will
% give a false positive is higher than any of the individual tests. We're
% doing something slightly different, which is asking whether 5 separate
% tests are consistent with an underlying phenomenon: that the person is
% "anxious."

% Q18: In the data provided in the .mat file, is there evidence that
% the tests results are correlated across patients?

% calculate correlations across patients
[rhoExpt,pValExpt] = corr(expt);
[rhoCtl,pValCtl] = corr(ctl);

figure('Name','Test Score Correlations');
subplot(2,1,1);
imagesc(rhoExpt,[0,1]);
colorbar; colormap('hot');
xlabel('Test #'); ylabel('Test #');
title('Patients');

subplot(2,1,2);
imagesc(rhoCtl,[0,1]);
colorbar; colormap('hot');
xlabel('Test #'); ylabel('Test #');
title('Controls');

% A18: Yes. In the experimental group, scores on all tests are correlated,
% with correlation coefficients ranging from 0.54 to 0.82, and all are
% statistically significant at p < 0.05. Similar for the control group,
% thought tests 1 and 5 do not appear significantly correlated:


% Q19: How should your findings about the correlations of the tests adjust
% the way you think about correcting for the fact that you're performing
% multiple comparisons?

% A19: They are at least partially independent tests, so the probability of
% observing at least one false positive across the whole set is higher than
% the probability of rejecting the null hypothesis for any single test. We
% should thus use a more stringent criterion for significance than p <
% 0.05, although not as stringent as if they were independent.

% Q#: The authors cite your work in the statistics section as a
% justification for not correcting for multiple comparisons. Based on their
% evidence, should they be able to conclude that there are two types of
% anxiety?

% A#: No. The two situations are quite different. Your paper asked "do the
% data support an underlying state of anxiety (as reflected by these 5
% measures)?" and so your hypothesis testing produced results which were
% correlated with each other and addressed a common question. These authors
% ask "Is there a detectable state of anxiety, as reflected by at least one
% test showing a significant difference?" and perform many hypothesis tests
% which produce independent conclusions. As a result, they should correct
% for false positives in those independent tests.

%% Simulate data with different levels of correlation:

% The data provided in 'anxiety_assays.mat' were made up. Data for each
% group was generated using the MATLAB function mvnrnd(mu,sigma,20), with a
% different 'mu' for the experimental and control groups. Specifically,
% muExpt = [1,1,1,1,1] and muCtl = [0,0,0,0,0]. They both used the same
% covariance matrix.
nSamp = 20;     % number per group
myCorr = 0.6;   % correlation among tests
muExpt = [1,1,1,1,1];   % normalized mean of experimental group
muCtrl = [0,0,0,0,0];   % normalized mean of control group
% covariance matrix:
mySigma = (ones(5,5) .* myCorr) + (eye(5) .* (1 - myCorr));

% Q20: What is the probability of at least one false positive?
% Key is to draw both samples from the SAME distribution.
nSim = 10000;
nFP = 0;
rng default
for k = 1:nSim
    R1 = mvnrnd(muCtrl,mySigma,nSamp);
    R2 = mvnrnd(muCtrl,mySigma,nSamp);
    H = ttest2(R1,R2);
    if any(H)
        nFP = nFP + 1;
    end
end
pFP = nFP / nSim;
% A20: 0.1694

% Q21: What is the probability of producing a false positive if the tests
% are completely uncorrelated (i.e. sigma is 0 everywhere except the
% main diagonal)?

% A21: There are two different ways to answer this question. The simplest
% is to just generate a new covariance matrix and repeat our simulation:
mySigma = eye(5);   % and repeat lines 119-130
% This gives 0.2279

% But it is easy to calculate exactly with independent tests, since the
% probability of one or more FPs is just 1 - all True Positives:
FWERuncorr = 1 - 0.95^nTests;
% This gives 0.2262

% Q22: What is the probability if all tests scores are perfectly correlated?

% A22: Again, two possible approaches. Simulation:
mySigma = ones(5,5);    % and repeat lines 119-130
% This gives 0.0504
% Or you can just realize that if the scores all are perfectly correlated,
% we have really only done one test, so FWER = alpha for single test.
% (Brian Healy gave an example in his multiple comparisons lecture of
% comparing heights between two groups, but coding the heights (single
% measurement) in five different ways: cm., in., ft., m., yds.