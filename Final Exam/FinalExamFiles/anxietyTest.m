% anxietyTest.m
%
% solution code for Thornquist test questions
% RTB wrote it, 27 Oct 2017

% Questions 15-22 will be related to the same data ("anxiety.mat" in the
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

% NOTE: originally called 'FWERwithCorrSim.m'

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

%% Do separate 2-sample t-tests for each test

% Calculate unadjusted p-values:
!!! Your code here

% Q15: Without making any post-hoc adjustments for multiple comparisons, on
% which tests is there a significant difference between the "experimental"
% and "control" treatments (at a significance level of p < 0.05)? (Check
% all that apply)

% Q16:  Using the standard Bonferroni correction are all of the tests still
% significant at p < 0.05?


% Q17: Should you be concerned about multiple comparisons here?


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
!!! Your code here

% Q19: How should your findings about the correlations of the tests adjust
% the way you think about correcting for the fact that you're performing
% multiple comparisons?

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

% You can see that every test is positively correlated with every other
% test. Use a simulation to estimate the probability that you get a false
% positive on one or more of the five tests (i.e. the 'Family-wise error
% rate' or FWER) using this same covariance matrix and an 'n' of 20
% patients in each group.

nSim = 10000;
nFP = 0;
rng default
for k = 1:nSim
    !!! Your code here
end

% Q21: What is the probability of producing a false positive if the tests
% are completely uncorrelated (i.e. sigma is 0 everywhere except the
% main diagonal)?
!!! Your code here

% Q22: What is the probability if all tests scores are perfectly correlated?
!!! Your code here
