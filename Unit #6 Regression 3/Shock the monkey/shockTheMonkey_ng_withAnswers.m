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

% Load the file into a table structure called 'ds'
fileName ='es5bRaw.xlsx';   % big effect
% fileName = 'js25aRaw.xlsx'; % small effect
ds = readtable(fileName);

% Each row is a single trial during which the monkey viewed a stochastic
% motion display whose signal strength is varied systematically (Coh) and
% whose direction is chosen according to the preferred direction of neurons
% at the stimulation site. At the end of each trial, the monkey chose one
% of two possible directions: preferred direction (PDchoice = 1) or the
% null direction (PDchoice = 0). Positive values of Coh indicate a stimulus
% moving in the preferred direction; negative values indicate null
% direction motion. On any given trial, microstimulation was applied to the
% electrode with a probability of 0.5 (Mstim = 1 on microstim trials; 0 on
% control trials).
%
% Columns:
% 1) Mstim: 0/1 for absence/presence of microstimulation
% 2) Coh: Signed strength of the visual stimulus (% coherent motion)
% 3) PDchoice: 0/1 for monkey choice in the null/pref direction
%
% Our scientific questions are "Did the microstimulation influence the
% monkey's perceptual decisions? If so, by how much?"

%% Plot the data

% TODO: We want a plot like that shown in the accompanying figure (see Q1
% of the Learning Catalytics module "Shock the Monkey"
% x-axis is the signed correlation value; y-axis is proportion of preferred
% decisions. Stim trials in red; no-stim trials in black

% define some constants for our column names
MSTIM=1; COH=2; PROP=3; CILO=4; CIHI=5;

% We first need to find all stimulus conditions:
allDataProp = sortrows(unique([ds.Mstim,ds.Coh],'rows'),[1,2]);
nConds = length(allDataProp);
% Add a column of zeros to allDataProp to hold the calcualted proportion;
% (plus two more for confidence intervals)
allDataProp = [allDataProp,zeros(nConds,3)];
myAlpha = 0.05;   % for error bars of the 95% CI

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
%
% ANSWER: length(ds.Mstim) = 560

% QUESTION (Q2): How many unique types of trial were there?
%
% ANSWER: nConds = 14

% QUESTION (Q3): How many repetitions of each trial type did the monkey
% perform?
%
% ANSWER: 560 / 14 = 40

% QUESTION (Q4): The error bars in the figure correspond to the 95% CI from
% the binomial distribution. Name two different ways in which we could make
% the error bars smaller?
%
% ANSWER:
% 1. Have the monkey perform more repetitions for each trial type.
% 2. Use a lower level of confidence for the CI.

%% Bonus: bootstrapping confidence intervals

% Suppose you didn't know that MATLAB had a handy-dandy function to
% generate confidence intervals. How might you get them? You should know
% this from the ASA/stroke exercise!

% Let's use an example a single trial type in which the monkey performed 40
% repetitions and chose the PD target on 16 of them:

% Answer using 'binofit':
[pHat,pCI] = binofit(16,40,myAlpha);
% pHat = 0.4000
% pCI = 0.2486    0.5667

% Answer using the bootstrap (percentile method)
% create a vector of 16 ones and 24 zeros:
x = [ones(16,1);zeros(40-16,1)];
% create an anonymous function that calculates the proportion:
prop = @(n) sum(n)/length(n);
% use 'bootci'
ci = bootci(10000,{prop,x},'alpha',myAlpha,'type','percentile');
% ci = 0.25, 0.55

%% Fit full model using glmfit

% TODO: Write down the full regression model, including an interaction term
% for microstimulation and signal strength.

% ln(P/1-P) = b0 + b1*stim + b2*corr + b3*stim*corr
[b,dev,stats] = glmfit([ds.Mstim,ds.Coh,ds.Mstim .* ds.Coh],ds.PDchoice,...
    'binomial','link','logit');

% QUESTION (Q5): Is the interaction term statitically significant at p < 0.05?
%
% ANSWER: stats.p(4) = 0.89
% Not significant.

% TODO: If your answer to the question is 'no', re-do the fit without the
% interaction term.
% The significance test for the interaction term is in stats.p(4)
if stats.p(4) > 0.05
    [b,dev,stats] = glmfit([ds.Mstim,ds.Coh],ds.PDchoice,...
        'binomial','link','logit');
end

% QUESTION (Q6): What is the scientific meaning of beta0 (i.e. the
% y-intercept) in our model?
%
% ANSWER: It is the log-odds of the monkey making a preferred-direction
% choice on trials where there is no microstimulation and the motion
% strength (coh) is 0.

% QUESTION (Q7): Look at your graph. What would a value of 0 for the beta0
% coefficient correspond to in terms of the probability of the monkey
% making a preferred decision choice on trials where there was no
% microstimulation and the motion strength was zero?
%
% ANSWER: Look at the regression equation! For trials without
% microstimulation with a motion strength of 0, the regression reduces to:
% log(P/1-P) = 0, or (P/1-P) = 1 and P = 0.5

% QUESTION (Q8): Did microstimulation affect the monkey's choices at a
% significance level of p < 0.05?
%
% ANSWER: Yes! From looking at the plots of the data and from the
% regression we can see that the effect is both large (i.e. perceptually
% significant) and highly statistically significant.

% QUESTION (Q9): What is the p-value for the model parameter capturing the
% effect of microstimulation on the monkey's choices?
%
% ANSWER: stats.p(2) = 5.0321e-20

%% Plot the model fits

% TODO: Plot the regression lines on top of the raw data. Make a separate
% line for the stim (red line) and no-stim (black line) predictions.
% Remember that our y-axis is in units of 'proportion preferred decisions'
% and NOT in log(P/1-P). HINT: You need to solve for 'P'

% Range for coherence:
coh = min(ds.Coh):0.01:max(ds.Coh);
const = ones(size(coh));
mStim = ones(size(coh));
pStim = 1 ./ (1 + exp(-(b(1).*const + b(2).*mStim + b(3).*coh)));
plot(coh,pStim,'r-');
pNoStim = 1 ./ (1 + exp(-(b(1).*const + b(2).*~mStim + b(3).*coh)));
plot(coh,pNoStim,'k-');

%% TODO: Determine the equivalent visual stimulus for microstimulation

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
% at which the animal is equally likely to report pref. or null. This is
% referred to in the psychophysical literature as the "Point of Subjective
% Equality." The key mathematical point is that this is where the lefthand
% side of our regression equation, log(p/1-p), is equal to 0. So we
% calculate the signal strength at PSE for the stim curve and subtract this
% from the signal strength at PSE for the ctrl curve. This will give us an
% equation in terms of the beta parameters in our model. Very elegant!
EVS1 = b(2) / b(3);
eqv = round(EVS1 * 10) / 10;
tStr = [fileName ': Equiv. Visual Stimulus = ',num2str(eqv), '% Coh'];
title(tStr);

% 2) Brute force. We have the regression equation, so we can input a very
% finely spaced set of coh values and find the one that gives us a value of
% 0.5. Since we're unlikely to get exactly 0.5, we would choose some narrow
% range straddling 0.5 and then take the average.
EVS2 = mean(coh(pNoStim < 0.505 & pNoStim > 0.495)) - ...
      mean(coh(pStim < 0.505 & pStim > 0.495));
  
% We can encourage the students to break this down into steps. First create
% a logical selection vector that has ones where the y-value (pNoStim or
% pStim) are within some narrow range around 0.5. Use these to find the
% corresponding x-values (in coh), then take the average of these values.
PSEnoStim = mean(coh(pNoStim < 0.505 & pNoStim > 0.495));
PSEstim = mean(coh(pStim < 0.505 & pStim > 0.495));

% QUESTION (Q10): How much signal would need to be added to the random dot
% display in order to match the effect of microstimulation on the monkey's
% choices? Give your answer as a positive percentage to 1 decimal place.

% QUESTION (Q11): Upload your final figure to Learning Catalytics.

% Plot the PSE lines:
ax = axis;
line([ax(1),ax(2)],[0.5,0.5],'Color','k','LineStyle','--');
line([PSEstim,PSEstim],[ax(3),0.5],'Color','r','LineStyle','--');
line([PSEnoStim,PSEnoStim],[ax(3),0.5],'Color','k','LineStyle','--');
p1 = plot(PSEstim,0,'r.','MarkerSize',15);
p2 = plot(PSEnoStim,0,'k.','MarkerSize',15);

% add a legend'
legend([h1,h2,p1,p2],{'Stim','No Stim','PSE_s','PSE_n_s'},'Location','Northwest');