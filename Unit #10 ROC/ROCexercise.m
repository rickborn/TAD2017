% ROCexercise.m
%
% Signal Detection Theory with spikes from a retinal ganglion cell
%
% RTB wrote it, 03 December 2017
%
% adapted from NB204_ROC_Demo.m by RTB on 14 April 2003

% Homework instructions:
%
% Collaboration and consultation of outside sources. The goal of this
% exercise is to give you practice in plotting and understanding ROC curves
% in MATLAB. The instructions and supplied code below will lead you through
% this in a straightforward and gentle manner. Therefore, we ask that you
% do not collaborate with any classmates on this exercise prior to class
% discussion. Neither should you seek out answers from published materials
% or online code repositories. This is a self-contained exercise, so you
% should just work through it on your own.
%
% What to do: Each section begins with a '%%'. Read through the commented
% code and follow the instructions provided. In some cases you will be
% asked to answer a question, clearly indicated by 'QUESTION'--please do so
% using comments (i.e. by starting each line with a '%'. In other cases,
% you will be asked to supply missing code, indicated by 'TODO'. Once you
% have supplied the required code, you can execute that section by
% mouse-clicking in that section (The block will turn yellow.) and then
% then simultaneously hitting the 'ctrl' and 'enter' keys (PC) or 'command'
% and 'enter' keys (Mac).
%
% The first section of code doesn't require you to do anything, but you
% still need to execute it in order to load the data and put up an
% initial plot.

%% Example from Barlow, Levick & Yoon 1971
% (fig. 3 from Parker & Newsome 1998)

% Barlow et al (1971) measured the responses of cat retinal ganglion cells
% to brief flashes of light. They obtained complete probability
% distributions for the occurrence of 0, 1, 2 . . . N spikes in a 200-ms
% counting window for stimulus intervals that had weak flashes of light (~5
% photons, 'signal') and for equivalent blank intervals (0 photons,
% 'noise').
%
% We want a metric of how well our neuron can distinguish signal from
% noise.

%% Load data and plot a histogram of responses to signal vs. noise

% The data file consists of spike counts measured during the 200-ms
% counting window. Column #1 is for the trials containing signal (~5
% photons); column #2 is for trials containing only noise (0 photons).
load BLYdata

% We want to have 1 bin for each possbile spike count
xBins = min(spikeCounts(:)):1:max(spikeCounts(:));
figure
subplot(3,1,1);
hist(spikeCounts,xBins);
xlabel('Response (spikes/trial)');
ylabel('# of trials');
legend('signal','noise');
title('Spike count histograms');

set(gcf,'Position',[456 29 837 779],'Name','Barlow, Levick & Yoon 1971');

% QUESTION:
% What are the average spike counts for signal and noise trials?
% ANSWER:
% mean(spikeCounts);      % 6.6354    3.9583

%% Performance metric #1: Joint probabilities

% 2AFC = Two-alternative forced choice. That is, on each trial, we have two
% possible responses ('flash YES' or 'NO flash'), and we must give one
% answer on each trial.
%
% Now consider a 2AFC game in which we are presented, on each “trial,” with
% two spike counts, one of which was drawn from the signal (blue)
% distribution and the other of which was drawn from the noise (yellow)
% distribution, and we are asked to decide which measurement corresponds to signal
% (i.e. came from the blue distribution).
% 
% Looking at the two distributions, the only reasonable decision strategy
% would be to assign the larger response to the signal (5 quanta) and the
% smaller response to the noise (0 quanta). If both responses are the same
% size, we have to guess (for instance, by tossing a coin).
% 
% How often will we get the correct answer? We can answer this precisely by
% looking at joint probabilities. We know from our data what the
% probability is of getting, for instance, a response of exactly 1 spike on
% trials where a dim flash (5 quanta of light) was presented: It happened
% in five trials out of 96, so the probability is 0.0521. If we have
% two measurements, we are interested in the joint probability: for
% instance, what is the probability of getting 4 spikes from the blue
% distribution AND 2 spikes from the yellow distribution? Since the draws
% are independent, we can compute the joint probabilities by multiplying
% the single probability distributions.

% Use hist to return the counts in each bin:
[countsPerBin,~] = hist(spikeCounts,xBins);

% TODO: Calculate the probabilities of each spike count for each condition.
% HINT: For each distribution, the probabilities must sum to 1.
pSignal = countsPerBin(:,1) ./ sum(countsPerBin(:,1));
pNoise =  countsPerBin(:,2) ./ sum(countsPerBin(:,2));

% QUESTION:
% What is the probability of getting 0 spikes on a 'signal' trial?
% ANSWER:
% pSignal(1);   % 0.0104

% QUESTION:
% What is the probability of getting 0 spikes on a 'noise' trial?
% ANSWER:
% pNoise(1);   % 0.0313

% TODO: Calculate the joint probability distribution and store it in 'J'
% HINT: Our possible spike counts can range from 0 to 17, so 'J' should be
% an 18 x 18 matrix, where bin i,j contains the probability of getting i-1
% spikes from the signal distribution AND j-1 spikes from the noise
% distribution. So, for example, J(1,1) should contain the probability of
% getting 0 spikes on a signal trial AND 0 spikes on a noise trial.
% NOTE: The reason for the weird offset (e.g. j-1 spikes) is that we are
% considering spike counts from 0 to 17, but MATLAB indices can't be less
% than 1. So an index of 1 corresponds to a spike count of 0 and an index
% of 18 corresponds to a spike count of 17.
J = pSignal * pNoise';

% QUESTION: If we take a random draw from each distribution, what is the
% probability of getting 0 spikes from the 'signal' distribution AND 0
% spikes from the 'noise' distribution?
% ANSWER:
% pSignal(1) * pNoise(1);   % 3.2552e-04
% J(1,1);   % 3.2552e-04

% Plot the join probability distribution
subplot(3,1,2);
imagesc(J);
xlabel('P(n spikes) from noise'); ylabel('P(n spikes) from signal');
title('Joint probability distribution');
c = colorbar;
c.Label.String = 'Probability';

%% Calculate the probability of being correct from the table:

% We can now read all of the joint probabilities directly from this table.
% For instance, we can look up the probability of getting 4 spikes from the
% blue distribution AND 2 spikes from the yellow (in this case, the
% probability is around 0.02).
% 
% In order to see whether our decision rule (assign the bigger number of
% spikes to the signal) is any good, we need to determine how likely it is
% to get more spikes with signal than with noise. That is, we are
% interested in the lower triangular part of the joint-probability table.
% 
% Summing over all of these will give us the total probability of getting a
% higher spike rate from the signal than the noise distribution.
% 
% In addition, if both samples yield the same spike count, we have to assign
% them to signal vs. noise at random, and if we do that, we will be
% correct, on average, 50% of the time. So, to the sum of the probabilities
% in the lower triangle, we add 0.5 times the sum of the probabilities along
% the main diagonal and thus get the overall probability of getting it
% right.

% TODO: Calculate the probability of being correct
% HINT: Look up the MATLAB documentation for 'tril' and 'diag'
pCorrectJP = sum(sum(tril(J,-1))) + 0.5*(sum(diag(J,0)));

% QUESTION: What is the value of 'pCorrectJP'?
% ANSWER: 0.7228

%% Performance method #2: Receiver Operating Characteristic

% We want to know how well an "ideal observer" ("ideal" because she has
% access to all of the statistics contained in the two histograms) could
% distinguish signal from noise. Conceptually this amounts to adopting
% different "criteria"—that is a number of spikes at or above which we
% would declare that a signal was present. For example, if we adopted an
% extremely strict criterion of 17 spikes we would correctly declare a
% "hit" (i.e. signal present) only 1 time out of the total of 96 blue
% trials, for a "hit rate" of 1/96, or 0.0104. We would never declare a
% signal to be present when it wasn't (a "false alarm"), because there are
% no yellow bars at or to the right of 17. This gives us two values, 0.0104
% and 0, which constitute one point on a plot of "hits" vs. "false alarms."
% We then repeat this procedure for all possible criteria from 17 to 0,
% each time sliding one to the left and calculating our "hits" from the
% signal distribution and our "false alarms" from the noise distribution.

maxVal = max(spikeCounts(:));       % maximum spike count
minVal = min(spikeCounts(:));       % minimum spike count
totTrials = length(spikeCounts);    % # of trials (n = 96 per condition)
allCriteria = [maxVal:-1:minVal]';  % all criteria that will be considered
pHit = zeros(size(allCriteria));    % variable to store hit rates
pFP = zeros(size(allCriteria));     % variable to store false positive rates

% TODO: For each spike-count criterion in 'allCriteria', calculate the
% probability of a hit ('pHit') and a false alarm ('pFA')
for k = 1:size(allCriteria,1)
    pHit(k) = sum(spikeCounts(:,1) >= allCriteria(k)) / totTrials;
    pFP(k) = sum(spikeCounts(:,2) >= allCriteria(k)) / totTrials;
end

%% Plot our ROC curve

subplot(3,1,3);
plot(pFP,pHit,'bs-');
hold on;
ylabel('Prob(Hit)');
xlabel('Prob(False Alarm)');
plot([0 1],[0 1],'k--');
title('ROC curve');

%% Calculate the area under the curve

% TODO: Calculate the area under our ROC curve
% HINT: Use 'trapz'
% If you want to understand how 'trapz' works, see:
% https://en.wikipedia.org/wiki/Trapezoidal_rule
pCorrectROC = trapz(pFP,pHit);

% QUESTION: What is the value of 'pCorrectROC'?
% ANSWER: 0.7228

% Print our two values on the final plot. How do they compare?
% If they are not identical, you did something wrong! Go back and check
% your ROC calculation. How did you count the # of hits and false alarms?
% Remember that we consider our criterion met if the value in the
% distribution is greater than OR equal to the criterion.
tStr1 = sprintf('sum JP = %.4f',pCorrectJP);
tStr2 = sprintf('auROC = %.4f',pCorrectROC);
text(0.6,0.3,tStr1);
text(0.6,0.2,tStr2);

% QUESTION: Does the ROC method make any assumptions about how our spike
% count values are distributed?
% ANSWER: No. ROC is a nonparametric method, so it is applicable to
% non-normally distributed data, which is often true of spike counts, which
% tend to follow a Poisson distribution.