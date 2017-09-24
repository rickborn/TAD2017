% mStimLogisticRegressionDemo.m
%
% This is a more tutorial-like version adapted from mstimfit.m and
% plotmstimfit.m
%
% Data a from Salzman et al. (1992) "Microstimulation in visual area MT:
% Effects on direction discrimination performance," J. Neurosci.
% 12(6):2331-2355.
%
% Note: In the paper, only proportions are given in the figures so I
% back-generated simulated raw data sets based on the reported # of trials
% for each condition. see: mstimdatagen.m
%
% RTB wrote it, 26 January 2017

% Concepts covered: 
% 1. Creating summary statistics (proportions) by matching rows in raw data 
% 2. Logistic regression with glmfit & fitglm
% 3. Indicator (dummy) variables for binomial covariates 
% 4. Model specification using Wilkinson notation
% 5. Model comparison: full(+interaction) vs. reduced
% 6. Scientific interpretation of model fits
% 7. Confidence intervals on model predictions

%% Read in the data from a Newsome lab microstimulation experiment
%fileName = 'es5bRaw.xlsx';
fileName = 'js25aRaw.xlsx';
ds = dataset('xlsfile',fileName);

% Each row is a single trial during which the animal viewed a stochastic
% motion display whose signal strength is varied systematically (Coh) and
% whose direction is chosen according to the preferred direction of neurons
% at the stimulation site. At the end of each trial, the monkey chose one
% of two possible directions: preferred direction (PDchoice = 1) or the
% null direction (PDchoice = 0). Positive values of Coh indicate a stimulus
% moving in the preferred direction; negative values indicate null
% direction motion. On any given trial, microstimulation was applied to the
% electrode with a probability of 0.5 (Mstim = 1 on microstim trials; 0 on
% control trials).
% Columns:
% 1) Mstim: 0/1 for absence/presence of microstimulation
% 2) Coh: Stimulus coherence
% 3) PDchoice: 0/1 for monkey choice in the null/pref direction
%
% To see these: fieldnames(ds)
%
% Our scientific question is "Did the microstimulation influence the
% monkey's perceptual decisions? If so, by how much?"
%% Plot the data
% We want a plot like those in fig. 4 of Salzman et al. 1992:
% x-axis is the signed correlation value; y-axis is proportion of preferred
% decisions.

% We first need to find all stimulus conditions:
allDataProp = sortrows(unique([ds.Mstim,ds.Coh],'rows'),[1,2]);
nConds = length(allDataProp);
% Add a column of zeros to allDataProp to hold the calcualted proportion
allDataProp = [allDataProp,zeros(nConds,1)];
for k = 1:nConds
    thisCond = find(ds.Mstim ==allDataProp(k,1) & ds.Coh == allDataProp(k,2));
    allDataProp(k,3) = sum(ds.PDchoice(thisCond)) / length(thisCond);
end

% Now plot: stim trials in red; no-stim trials in black
stimIdx = logical(allDataProp(:,1));
figure
plot(allDataProp(stimIdx,2),allDataProp(stimIdx,3),'ro','MarkerFaceColor','r');
hold on;
plot(allDataProp(~stimIdx,2),allDataProp(~stimIdx,3),'ko');
xlabel('Motion Strength (%coh)'); ylabel('Proportion Preferred Decisions');
ax = axis;
axis([ax(1),ax(2),0,1]);
legend('Stim','No Stim','Location','Northwest');
title(fileName);

%% Fit full model using glmfit
% ln(P/1-P) = b0 + b1*stim + b2*corr + b3*stim*corr
[b,dev,stats] = glmfit([ds.Mstim,ds.Coh,ds.Mstim .* ds.Coh],ds.PDchoice,'binomial','link','logit');

% The significance test for the interaction term is in stats.p(4)
if stats.p(4) > 0.05
    [b,dev,stats] = glmfit([ds.Mstim,ds.Coh],ds.PDchoice,'binomial','link','logit');
end

%% Fit full model using fitglm

% Search MATLAB documentation for 'Wilkinson Notation'
modelspec = 'PDchoice ~ 1 + Mstim*Coh'; % A*B is same as A + B + A:B
mdlFull = fitglm(ds,modelspec,'Distribution','binomial');

% Compare 'b' with mdlFull.Coefficients
if mdlFull.Coefficients{4,'pValue'} > 0.05
    modelspec = 'PDchoice ~ 1 + Mstim + Coh';
    mdlReduced = fitglm(ds,modelspec,'Distribution','binomial');
end

%% Plot the model fits

coh = min(ds.Coh):0.01:max(ds.Coh);
const = ones(size(coh));
mStim = ones(size(coh));
pNoStim = 1 ./ (1 + exp(-(b(1).*const + b(2).*~mStim + b(3).*coh)));
plot(coh,pNoStim,'k-');
pStim = 1 ./ (1 + exp(-(b(1).*const + b(2).*mStim + b(3).*coh)));
plot(coh,pStim,'r-');

%% Determine the equivalent visual stimulus for microstimulation

% Two possible methods:
% 1) Using the regression equation, solve for the signal strength at which
% the animal is equally likely to report pref. or null. This is referred to
% in the psychophysical literature as the "Point of Subjective Equality."
% The key mathematical point is that this is where the lefthand side of our
% regression equation, log(p/1-p), is equal to 0. So we calculate the
% signal strength at PSE for the stim curve and subtract this from the
% signal strength at PSE for the ctrl curve. This ends up to be b(2)/b(3)
equivVisualStimulus = b(2) / b(3);
eqv = round(equivVisualStimulus * 10) / 10;
tStr = [fileName ': Equiv. Visual Stimulus = ',num2str(eqv), '%coh'];
title(tStr);

% 2) Brute force approach. We have the regression equation, so we can input
% a very finely spaced set of coh values and find the one that gives us a
% value of 0.5. Since we're unlikely to get exactly 0.5, we would choose some
% narrow range straddling 0.5 and then take the average. 
EVS = mean(coh(pNoStim < 0.505 & pNoStim > 0.495)) - ...
      mean(coh(pStim < 0.505 & pStim > 0.495));

% We can encourage the students to break this down into steps. First create
% a logical selection vector that has ones where the y-value (pNoStim or
% pStim) are within some narrow range around 0.5. Use these to find the
% corresponding x-values (in coh), then take the average of these values.
PSEnoStim = mean(coh(pNoStim < 0.505 & pNoStim > 0.495));
PSEstim = mean(coh(pStim < 0.505 & pStim > 0.495));

% Draw lines for the respective PSEs
ax = axis;
line([ax(1),PSEstim],[0.5,0.5],'Color','r','LineStyle','-');
line([PSEstim,PSEstim],[ax(3),0.5],'Color','r','LineStyle','--');
line([ax(1),PSEnoStim],[0.5,0.5],'Color','k','LineStyle','--');
line([PSEnoStim,PSEnoStim],[ax(3),0.5],'Color','k','LineStyle','--');

 %% Using the 'predict' function, we can get CIs on our estimates
 
 % Re-plot raw data
figure
plot(allDataProp(stimIdx,2),allDataProp(stimIdx,3),'ro','MarkerFaceColor','r');
hold on;
plot(allDataProp(~stimIdx,2),allDataProp(~stimIdx,3),'ko');
xlabel('Correlation (%)'); ylabel('Proportion Preferred Decisions');
ax = axis;
axis([ax(1),ax(2),0,1]);
title(tStr);

% Get predictions and CIs for stim trials
coh = min(ds.Coh):max(ds.Coh);
mStim = ones(size(coh));
[yStim,yStimCI] = predict(mdlReduced,[mStim',coh']);
% Plot 'em
plot(coh',yStim,'r-','LineWidth',2);
plot(coh',yStimCI(:,1),'r--');
plot(coh',yStimCI(:,2),'r--');

% Get predictions and CIs for stim trials
[yNoStim,yNoStimCI] = predict(mdlReduced,[~mStim',coh']);
% Plot 'em
plot(coh',yNoStim,'k-','LineWidth',2);
plot(coh',yNoStimCI(:,1),'k--');
plot(coh',yNoStimCI(:,2),'k--');