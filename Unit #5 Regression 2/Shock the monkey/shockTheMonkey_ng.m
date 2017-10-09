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

% TODO: Load the file 'js25aRaw.xlsx' in a dataset structure called 'ds'
ds = ;

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

!!! Your code here
plot(..., 'ro','MarkerFaceColor','r');  % trials WITH microstimulation
hold on
plot(..., 'ko');                        % trials WITHOUT micrsotimulation
xlabel('Motion Strength (%coh)'); ylabel('Proportion Preferred Decisions');

%% Fit full model using glmfit

% TODO: Write down the full regression model, including an interaction term
% for microstimulation and signal strength.

% TODO: Use 'glmfit' to perform the regression
[b,dev,stats] = glmfit();

% QUESTION (Q1): Is the interaction term statitically significant at p < 0.05?

% TODO: If your answer to the question is 'no', re-do the fit without the
% interaction term.

%% Plot the model fits

% TODO: Plot the regression lines on top of the raw data. Make a separate
% line for the stim (red line) and no-stim (black line) predictions.
% Remember that our y-axis is in units of 'proportion preferred decisions.'

!!! Your code here

% QUESTION (Q2): Did microstimulation affect the monkey's choices at a
% significance level of p < 0.05?

% QUESTION (Q3): What is the p-value for the model parameter capturing the
% effect of microstimulation on the monkey's choices?


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

% 2) Brute force. We have the regression equation, so we can input a very
% finely spaced set of coh values and find the one that gives us a value of
% 0.5. Since we're unlikely to get exactly 0.5, we would choose some narrow
% range straddling 0.5 and then take the average.

!!! Your code here

% QUESTION (Q4): How much signal would need to be added to the random dot
% display in order to match the effect of microstimulation on the monkey's
% choices? Give your answer as a positive percentage to 1 decimal place.

% QUESTION (Q5): Upload your final figure to Learning Catalytics.
