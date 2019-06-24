% shockTheMonkey_ng.m
%
% Data a from Salzman et al. (1992) "Microstimulation in visual area MT:
% Effects on direction discrimination performance," J. Neurosci.
% 12(6):2331-2355.
%
% RTB wrote it, 26 January 2017
% RTB adapted for self-test, 09 October 2017, rainy day

% What to do: Login to Learning Catalytics (LC) and join the session for
% the module entitled "Shock the Monkey". You will answer a series of
% questions based on the guided programming below. Each section begins with
% a '%%'. Read through the comments and follow the instructions provided.
% In some cases you will be asked to answer a question, clearly indicated
% by 'QUESTION' and a corresponding 'Q#' that directs you to answer the
% relevant question in LC. In other cases, you be asked to supply missing
% code, indicated by 'TODO'. Once you have supplied the required code, you
% can execute that section by mouse-clicking in that section (The block
% will turn yellow.) and then simultaneously hitting the 'ctrl' and 'enter'
% keys (PC) or 'command' and 'enter' keys (Mac).

%% Read in the data from a Newsome lab microstimulation experiment

% Load a data file into a table structure called 'ds'
fileName ='es5bRaw.xlsx';
ds = readtable(fileName);

% Each row is data from a single trial during which the monkey viewed a
% stochastic motion display whose signal strength was varied systematically
% (Coh) and whose direction was chosen according to the preferred direction
% of neurons at the stimulation site. At the end of each trial, the monkey
% chose one of two possible directions: preferred direction (PDchoice = 1)
% or the null direction (PDchoice = 0). Positive values of Coh indicate a
% stimulus moving in the preferred direction; negative values indicate null
% direction motion. On any given trial, microstimulation was applied to the
% electrode with a probability of 0.5 (Mstim = 1 on microstim trials; 0 on
% control trials).
%
% Columns:
% 1) Mstim: 0/1 for absence/presence of microstimulation
% 2) Coh: Signed strength of the visual stimulus (% coherent motion)
%    positive values indicate motion in the neuron's 'preferred' direction;
%    negative values correspond to the opposite, or 'null', direction
% 3) PDchoice: 0/1 for monkey choice in the null/preferred direction
%
% Our scientific questions are "Did the microstimulation influence the
% monkey's perceptual decisions? If so, by how much?"

%% Plot the data

% TODO: We want a plot in which the x-axis is the signed correlation value;
% y-axis is proportion of preferred decisions. Stim trials in red; no-stim
% trials in black. NOTE: Because this step is rather time-consuming, I've
% provided the code here, so you don't need to do anything except execute
% it to get the figure. But please read through the code and make sure you
% understand how it works. And there are LC questions to be answered at the
% bottom of the cell.

% define some constants for our column names
MSTIM=1; COH=2; PROP=3; CILO=4; CIHI=5;

% We first need to find all of the different stimulus conditions:
allDataProp = sortrows(unique([ds.Mstim,ds.Coh],'rows'),[1,2]);
nConds = length(allDataProp);

% Add a column of zeros to allDataProp to hold the calcualted proportion;
% (plus two more for confidence intervals)
allDataProp = [allDataProp,zeros(nConds,3)];
myAlpha = 0.05;   % for error bars of the 95% CI

% For each unique condition, we need to find all of the corresponding rows
% in our original raw data file and tally the proportion of those trials on
% which the monkey chose the preferred direction choice target:
for k = 1:nConds
    thisCond = find(ds.Mstim ==allDataProp(k,MSTIM) & ds.Coh == allDataProp(k,COH));
    [pHat,pCI] = binofit(sum(ds.PDchoice(thisCond)),length(thisCond),myAlpha);
    allDataProp(k,PROP) = pHat;
    allDataProp(k,[CILO,CIHI]) = pCI;
end
% logical for indexing:
stimIdx = logical(allDataProp(:,MSTIM));

% NOTE: The 'errorbar' function plots error bars that are L(i) + U(i) long.
% That is, it doesn't treat our CI as an interval, but rather as a distance
% from the mean to the end of each error bar. So we need to subtract each
% from the mean:
allDataProp(:,CILO) = allDataProp(:,PROP) - allDataProp(:,CILO); % lower error bar
allDataProp(:,CIHI) = allDataProp(:,CIHI) - allDataProp(:,PROP); % upper error bar

% Now plot: stim trials in red; no-stim trials in black
figure
h1=errorbar(allDataProp(stimIdx,COH),allDataProp(stimIdx,PROP),...
    allDataProp(stimIdx,CILO),allDataProp(stimIdx,CIHI),'ro');
hold on
h2=errorbar(allDataProp(~stimIdx,COH),allDataProp(~stimIdx,PROP),...
    allDataProp(~stimIdx,CILO),allDataProp(~stimIdx,CIHI),'ko');

xlabel('Motion strength (%coh)'); ylabel('Proportion preferred decisions');
ax = axis;
axis([ax(1),ax(2),0,1]);
legend([h1,h2],'Stim','No Stim','Location','Northwest');
title(fileName);
set(gca,'FontSize',12); % make text larger

% QUESTION (Q1): How many trials did the monkey perform in this experiment?

% QUESTION (Q2): How many unique types of trial were there?

% QUESTION (Q3): How many repetitions of each trial type did the monkey
% perform?

% QUESTION (Q4): The error bars in the figure correspond to the 95% CI from
% the binomial distribution. Name two different ways in which we could make
% the error bars smaller?

%% Fit full model using glmfit

% TODO: Write down the full regression model, including an interaction term
% for microstimulation and signal strength. Remember that you want to
% pass the RAW data as your arguments to 'glmfit'

% TODO: Use 'glmfit' to perform the regression
[b,dev,stats] = glmfit();

% QUESTION (Q5): Is the interaction term statitically significant at p < 0.05?

% TODO: If your answer to the question is 'no', re-do the fit without the
% interaction term.
%
% The significance test for the interaction term is in stats.p(4)
if stats.p(4) > 0.05
    [b,dev,stats] = glmfit();
end

% QUESTION (Q6): What is the scientific meaning of beta0 (i.e. the
% y-intercept) in our model?

% QUESTION (Q7): Look at your graph. What would a value of 0 for the beta0
% coefficient correspond to in terms of the probability of the monkey
% making a preferred decision choice on trials where there was no
% microstimulation and the motion strength was zero?

% QUESTION (Q8): Did microstimulation affect the monkey's choices at a
% significance level of p < 0.05?

% QUESTION (Q9): What is the p-value for the model parameter capturing the
% effect of microstimulation on the monkey's choices?

%% Plot the model fits

% TODO: Plot the regression lines on top of the raw data. Make a separate
% line for the stim (red line) and no-stim (black line) predictions.
% Remember that our y-axis is in units of 'proportion preferred decisions'
% and NOT in log(P/1-P). HINT: You need to solve for 'P':
% P = 1 / 1 + exp[-(b0 + b1*stim + b2*corr)]
!!! Your code here

%% TODO: Determine the "equivalent visual stimulus" for microstimulation

% We want to know not just "if" microstimulation had an effect; we would
% also like to estimate the magnitude of the effect. Look at the plot you
% made. You can see that there is a certain signal strength at which the
% monkey is equally likely to choose the preferred vs. the non-preferred
% direction of motion. To do this, we would draw a horizontal line from the
% y-axis value of 0.5 to our curve, then drop a vertical line to the x-axis
% to get a value of the motion strength. We call this value the "Point of
% Subjective Equality" (PSE), because it represents the visual stimulus for
% which the monkey was indifferent to the two choices. We can calculate the
% effect of microstimulation in units of the visual stimulus (%coh) by
% subtracting the PSE during microstimulation trials from that during
% control trials.

% Two possible approaches:

% 1) Algebra. Using the regression equation, solve for the signal strength
% at which the animal is equally likely to report preferred or null. This
% is referred to in the psychophysical literature as the "Point of
% Subjective Equality" or "PSE". We calculate the signal strength at PSE
% for the stim curve and subtract this from the signal strength at PSE for
% the ctrl curve. This will give us an equation in terms of the beta
% parameters in our model. Very elegant!

% 2) Brute force. We have the regression equation, so we can input a very
% finely spaced set of coh values and find the one that gives us a value of
% 0.5. Since we're unlikely to get exactly 0.5, we would choose some narrow
% range straddling 0.5 and then take the average.

!!! Your code here

% QUESTION (Q10): How much signal would need to be added to the random dot
% display in order to match the effect of microstimulation on the monkey's
% choices? Give your answer as a positive percentage to 1 decimal place.

% QUESTION (Q11): Upload your final figure to Learning Catalytics.

%% Bonus questions

% 1) The GLM is a beutiful framework, because we get so much "for free,"
% such as standard errors on our betas as well as significance tests for
% whether or not they are different from 0. But let's suppose this wasn't
% the case. How would you design your own procedure to explicitly test our
% primary null hypothesis, which is that microstimulation has no effect on
% the monkey's choices?

% 2) We have a nice point estimate an "effect size" (or "equivalent visual
% stimulus") of 14.4% coh, but we would also like to know its precision.
% Design a procedure to determine the standard error for the effect size.
