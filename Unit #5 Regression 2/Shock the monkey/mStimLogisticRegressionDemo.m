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
fileName = 'es5bRaw.xlsx'; % big effect
%fileName = 'js25aRaw.xlsx'; % small effect; used in class
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
% Add a column of zeros to allDataProp to hold the calcualted proportion;
% (plus two more for confidence intervals)
allDataProp = [allDataProp,zeros(nConds,3)];
% myErr = 0.32;   % for error bars of the 68% CI (= SEM)
myErr = 0.05;   % for error bars of the 95% CI
for k = 1:nConds
    thisCond = find(ds.Mstim ==allDataProp(k,1) & ds.Coh == allDataProp(k,2));
    % allDataProp(k,3) = sum(ds.PDchoice(thisCond)) / length(thisCond);
    [pHat,pCI] = binofit(sum(ds.PDchoice(thisCond)),length(thisCond),myErr);
    allDataProp(k,3) = pHat;
    allDataProp(k,[4,5]) = pCI;
end

% Now plot: stim trials in red; no-stim trials in black
stimIdx = logical(allDataProp(:,1));
figure
plot(allDataProp(stimIdx,2),allDataProp(stimIdx,3),'ro','MarkerFaceColor','r');
hold on;
plot(allDataProp(~stimIdx,2),allDataProp(~stimIdx,3),'ko','MarkerFaceColor','k');
xlabel('Motion strength (%coh)'); ylabel('Proportion preferred decisions');
ax = axis;
axis([ax(1),ax(2),0,1]);
legend('Stim','No Stim','Location','Northwest');
title(fileName);

%% Add error bars

% Many students seem to think that you can't put error bars on a
% proportion: you divide one number by another, so how can you calculate a
% standard error? Well, you don't need to since we know our underlying
% probability model is binomial. We can use this to calculate confidence
% intervals directly. The intuition is that we are much more confident of a
% proportion of 0.8 for 80 of 100 than for 8 of 10.

% To do this, we use 'binofit', but it is much easier to incorporate this
% into the 'for' loop in the previous section were we calculate the sum.

% NOTE: The 'errorbar' function plots error bars that are L(i) + U(i) long.
% That is, it doesn't treat our CI as an interval, but rather as a distance
% from the mean to the end of each error bar. So we need to subtract each
% from the mean:

allDataProp(:,4) = allDataProp(:,3) - allDataProp(:,4); % lower error bar
allDataProp(:,5) = allDataProp(:,5) - allDataProp(:,3); % upper error bar

errorbar(allDataProp(stimIdx,2),allDataProp(stimIdx,3),...
    allDataProp(stimIdx,4),allDataProp(stimIdx,5),'r.');

errorbar(allDataProp(~stimIdx,2),allDataProp(~stimIdx,3),...
    allDataProp(~stimIdx,4),allDataProp(~stimIdx,5),'ko');

%% Bootstrap error bars: one example

% answer using binofit:
[pHat,pCI] = binofit(16,40);
% pHat = 0.4000
% pCI = 0.2486    0.5667

% bootstrap method
% create a vector of 16 ones and 24 zeros:
x = [ones(16,1);zeros(40-16,1)];
% create an anonymous function that calculates the proportion:
prop = @(n) sum(n)/length(n);
% use bootci
ci = bootci(10000,{prop,x},'type','percentile');
% ci = 0.25, 0.55

% pretty good, but at what cost?
tic;[pHat,pCI] = binofit(16,40);t1=toc;
tic;ci = bootci(10000,{prop,x},'type','percentile');t2=toc;
t2/t1

%% Plot the raw data, too?

nTrials = length(ds);
mstimIndex = logical(ds.Mstim);

% Will need to jitter data to see all of the points
jitterX = (rand(nTrials,1) * 3) - 1.5;
jitterY = (rand(nTrials,1) * 0.2) - 0.1;

plot(ds.Coh(mstimIndex)+jitterX(mstimIndex),ds.PDchoice(mstimIndex)+jitterY(mstimIndex),'r.');
plot(ds.Coh(~mstimIndex)+jitterX(~mstimIndex),ds.PDchoice(~mstimIndex)+jitterY(~mstimIndex),'k.');

axis([ax(1),ax(2),-0.1,1.1]);

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
pFlag = 0;
if pFlag
    ax = axis;
    line([ax(1),PSEstim],[0.5,0.5],'Color','r','LineStyle','-');
    line([PSEstim,PSEstim],[ax(3),0.5],'Color','r','LineStyle','--');
    line([ax(1),PSEnoStim],[0.5,0.5],'Color','k','LineStyle','--');
    line([PSEnoStim,PSEnoStim],[ax(3),0.5],'Color','k','LineStyle','--');
end

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