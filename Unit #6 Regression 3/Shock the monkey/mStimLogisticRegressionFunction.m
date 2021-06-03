function [b,stats,eqv] = mStimLogisticRegressionFunction(fileName,pFlag)

% mStimLogisticRegressionFunction: function version of mStimLogisticRegressionDemo.m
%
% [b,stats,eqv] = mStimLogisticRegressionFunction(fileName,pFlag)
%
% Inputs:
% - fileName, default 'js25aRaw.xlsx'
% - pFlag, = 1 to plot the data and the fits (default, no plots)
%
% Outputs:
% - b, the beta coefficients for the fit
%   b(1) is the intercept, modeling choice bias
%   b(2) models the effect of microstimulation
%   b(3) models the effect of the visual stimulus
%   b(4) models the interaction between microstim and the stimulus
% - stats, various statistics including stats.p for p-values
% - eqv, the equivalent visual stimulus
%
% RTB adapted this on Oct. 6, 2017, post Chenghua's tenure party

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

% Read in the data from a Newsome lab microstimulation experiment
%fileName = 'es5bRaw.xlsx';
% NOTE: For an example of a significant interaction term, use js92dRaw.xlsx
if nargin < 2, pFlag = 0; end
if nargin < 1, fileName = 'js25aRaw.xlsx'; end
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

if pFlag
    % Plot the data
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
end

% Fit model using glmfit
% ln(P/1-P) = b0 + b1*stim + b2*corr + b3*stim*corr
[b,~,stats] = glmfit([ds.Mstim,ds.Coh,ds.Mstim .* ds.Coh],ds.PDchoice,'binomial','link','logit');

% The significance test for the interaction term is in stats.p(4)
if stats.p(4) > 0.05
    [b,~,stats] = glmfit([ds.Mstim,ds.Coh],ds.PDchoice,'binomial','link','logit');
end

% Fit model using fitglm
% Search MATLAB documentation for 'Wilkinson Notation'
modelspec = 'PDchoice ~ 1 + Mstim*Coh'; % A*B is same as A + B + A:B
mdl = fitglm(ds,modelspec,'Distribution','binomial');

% Compare 'b' with mdl.Coefficients
if mdl.Coefficients{4,'pValue'} > 0.05
    modelspec = 'PDchoice ~ 1 + Mstim + Coh';
    mdl = fitglm(ds,modelspec,'Distribution','binomial');
end

if pFlag
    % Plot the model fits
    coh = min(ds.Coh):0.01:max(ds.Coh);
    const = ones(size(coh));
    mStim = ones(size(coh));
    pNoStim = 1 ./ (1 + exp(-(b(1).*const + b(2).*~mStim + b(3).*coh)));
    plot(coh,pNoStim,'k-');
    if length(b) == 4
        pStim = 1 ./ (1 + exp(-(b(1).*const + b(2).*mStim + b(3).*coh + b(4).*coh.*mStim)));
    else
        pStim = 1 ./ (1 + exp(-(b(1).*const + b(2).*mStim + b(3).*coh)));
    end
    plot(coh,pStim,'r-');
end

% Determine the equivalent visual stimulus for microstimulation
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

if pFlag
    tStr = [fileName ': Equiv. Visual Stimulus = ',num2str(eqv), '%coh'];
    title(tStr);
       
    % Using the 'predict' function, we can get CIs on our estimates
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
    [yStim,yStimCI] = predict(mdl,[mStim',coh']);
    % for standard error, use the 67% CI instead of the 95%
    %[yStim,yStimCI] = predict(mdl,[mStim',coh'],'alpha',0.33);
    % Plot 'em
    plot(coh',yStim,'r-','LineWidth',2);
    plot(coh',yStimCI(:,1),'r--');
    plot(coh',yStimCI(:,2),'r--');
    
    % Get predictions and CIs for stim trials
    [yNoStim,yNoStimCI] = predict(mdl,[~mStim',coh']);
    % Plot 'em
    plot(coh',yNoStim,'k-','LineWidth',2);
    plot(coh',yNoStimCI(:,1),'k--');
    plot(coh',yNoStimCI(:,2),'k--');
end