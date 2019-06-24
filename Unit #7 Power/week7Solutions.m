% week7Solutions.m
%
% Answers to the questions posed in the Learning Catalytics Module for
% week #7.
%
% RTB created this file on 19 September 2018

%% QUESTION (Q1)

% You have created a genetic model of Amyotrophic Lateral Sclerosis (ALS,
% a.k.a. "Lou Gehrig's Disease"), and you want to test for early signs of
% motor neuron degeneration. You decide to test 1-month old mice on the
% "rotarod": each mouse is placed on a rotating rod and the time it takes
% for the animal to fall off is recorded. You initially run two groups of
% 20 mice each, comparing homozygotes (Dz) to heterozygous littermates
% (Ctrl). A 2-sample t-test fails to attain significance. But it is not too
% far off, so you and your PI decide on the following strategy: Since the
% mice are abundant and the test is relatively easy to perform, you will
% add 10 new animals to each group and then re-check your statistical test.
% You will continue to do this until you get a significant result (p <
% 0.05!) or until you've reached 100 animals per group, in which case you
% will determine the experiment to have failed.
% 
% What is your true false positive rate? Give your answer in % to the
% nearest whole number.

rng default;

nSims = 10000;
myAlpha = 0.05;
nMax = 100;
nAddObs = 10;
nInit = 20;

% The logic for this simulation is that we generate all nMax pairs of
% observations at once. Then, in our inner 'for' loop where we run the
% t-tests, we just extend the the length of our index into 'allSims'.
FP = zeros(nSims,1);
allNdx = [nInit:nAddObs:nMax];
for jSim = 1:nSims
    allSims = randn(nMax,2);
    for k = 1:length(allNdx)
        if ttest2(allSims([1:allNdx(k)],1), allSims([1:allNdx(k)],2),'Alpha',myAlpha)
            FP(jSim) = 1;
            break
        end
    end
end
FPrate = (sum(FP) / nSims) * 100;
% ANSWER: 16.73, which rounds to 17

%% QUESTION (Q2)

% Chastened by Uri Simonsohn, you and your PI decide to play it straight.
% Before you do your experiment, you do some research and decide that the
% smallest effect size, measured as d-prime, that would be worth publishing
% is d' = 0.57. How many animals should you run in each group in order to
% have an 80% probability of detecting an effect of this size?

desiredPower = 0.80;
dPrimeSim = 0.57;
% In this case, we can use the formula:
n = sampsizepwr('t2',[1,1],1+dPrimeSim,desiredPower);
% ANSWER: 50

% But we could also get an answer via simulation:
% One strategy would be to first get a rough answer using coarse sampling,
% then zoom in around it.
coarseFlag = 0;
if coarseFlag
    nMin = 5; nMax = 105; sampInt = 10; nSims = 10000;
else
    nMin = 45; nMax = 55; sampInt = 1; nSims = 100000;
end
allN = [nMin:sampInt:nMax]';
allPowSim = zeros(length(allN),1);
rng default
for k = 1:length(allN)
    % Generate our samples:
    allSamples = randn(allN(k),2,nSims);
    % To simulate a real d-prime, we just add our dPrime to one of the samples.
    allSamples(:,2,:) = allSamples(:,2,:) + dPrimeSim;
    
    % Now run our 2-sample t-test:
    [h,~] = ttest2(squeeze(allSamples(:,2,:)),squeeze(allSamples(:,1,:)));
    % Calculate the power:
    allPowSim(k) = sum(h) / nSims;
end
nReqd = min(allN(allPowSim >= desiredPower));

% Plot the relationship between n and power
figure
plot(allN,allPowSim,'bo-');
hold on
xlabel('N');
ylabel('Power');
ax = axis;
h1 = line([ax(1),ax(2)],[desiredPower,desiredPower]);
set(h1,'Color','k','LineStyle','--')
title(['Simulated effect size, dPrime = ' num2str(dPrimeSim)]);
grid on

%% QUESTION (Q3)

% For the same experiment, how big would d-prime have to be for you to have
% the same power (i.e. an 80% chance of detecting the effect) with an n of
% only 20 in each group? Give your answer to 2 decimal places.

% In this case, we can do it by brute force by getting the power for a
% large range of possible values of d-prime, then finding the one that
% gives us the desired power.

dPrime = 0.50:0.001:1;
pwrout = sampsizepwr('t2',[0,1],dPrime,[],20);
dPrimeCrit = mean(dPrime(pwrout > 0.79 & pwrout < 0.81));

% ANSWER: 0.9090

%% QUESTION (Q4) Meta-analysis 1

% Suppose that a large number of labs have done the same experiment: 50
% animals per group with a "true" effect size of d' = 0.57. Further suppose
% that only the labs that got a statistically significant effect are going
% to get their paper published. If you now survey this literature and
% collect all of the reported p-values, what percentage of these should be
% p <= 0.01? Give your answer to the nearest whole number.

% Constants for now; would be passed as arguments to a function:
n = 50;
nSims = 100000;
dPrimeSim = 0.57;

% Generate our samples:
allSamples = randn(n,2,nSims);
% To simulate a real d-prime, we just add our dPrime to one of the samples.
allSamples(:,2,:) = allSamples(:,2,:) + dPrimeSim;

% Now run our t-test:
[h,pVals] = ttest2(squeeze(allSamples(:,2,:)),squeeze(allSamples(:,1,:)));

% Find the proportion of published p-values less-than-or-equal to 0.01:
perCentAtCrit = round((sum(pVals <= 0.01) / sum(h))*100);

% ANSWER: 73%

%% QUESTION (Q5) Meta-analysis 2

% What would your answer to question #4 be if there were really no effect?
% (i.e. d' = 0)

% You shouldn't need to run a simulation to get this one. Under H0, we
% would expect that 1% of all p-values would be less than or equal to 0.01.
% But we are only considering published results, so we would have already
% selected for p-values < 0.05, which, under H0, should be 5%. So we would
% expect 1%/5% = 0.20, or 20%.

% But you can also go ahead and run the simulation above for Q4 but setting
% 'dPrimeSim' to 0 and see that you get the same answer.

%% QUESTION (Q6): Publication bias and effect size

% As for Q4, suppose that a large number of labs have done the same
% experiment but this time with just 5 animals per group and a "true"
% effect size of d' = 0.5. Further suppose that only the labs that got a
% statistically significant effect are going to get their paper published.
% What will be the median *published* effect size? For this analysis,
% determine your fitness for publication using a 1-tailed (i.e. right tail)
% 2-sample t-test and a criterion for publication of p < 0.05.

% Constants for now; would be passed as arguments to a function:
n = 5;
nSims = 100000;
dPrimeSim = 0.5;

% Generate our samples:
allSamples = randn(n,2,nSims);
% To simulate a real d-prime, we just add our dPrime to one of the samples.
allSamples(:,2,:) = allSamples(:,2,:) + dPrimeSim;

% Now calculate our actual d-primes (spread due to sampling error)
allDprime = squeeze(diff(mean(allSamples)));    % SD = 1

% and then run our t-test to see who gets published:
[h,pVals] = ttest2(squeeze(allSamples(:,2,:)),squeeze(allSamples(:,1,:)),...
    'Tail','Right');

% median published effect size:
medPubDprime = median(allDprime(logical(h')));
% Ans.: 1.31
% The mean is 1.34

%% QUESTION (Q7): Publication bias, effect size, and power

% Based on your simulation for Q6, what is the statistical power for n=5
% and an effect size of dPrime = 0.5?

% In the above simulation, we generated our simulated populations with two
% different means, so we know that ALL samples are truly different. So
% power is very easy to calculate:
powerSim = sum(h) / nSims;
% Ans.: 0.1777

% NOTE: If we re-run the above simulation using a 2-tailed test, we get a
% simulated power value of 0.1075.
pwrout = sampsizepwr('t2',[0,1],dPrimeSim,[],5);
% Ans. 0.1077

% We could also compute the power using 'sampsizepwr' and specifying a
% one-tailed test:
pwrout = sampsizepwr('t2',[0,1],dPrimeSim,[],5,'Tail','right');
% Ans. 0.1788