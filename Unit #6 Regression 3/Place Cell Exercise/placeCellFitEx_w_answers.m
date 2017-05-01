% placeCellFitEx.m
% 
% Chapter 9, of "Case studies in neural data analysis" by Kramer & Eden 2016
%
% modified by RTB 9 Dec. 2016

% Concepts covered:
% 1. Working with spike data: times to indices
% 2. Occupancy normalized histogram for place fields
% 3. Using glmfit to form a Poisson point-process model
% 4. Model selection through analysis of residuals
% 5. Model comparison through measures of goodness of fit: AIC,
%       Chi-square, parameter CIs, Kolmogorov-Smirnov 
%
% Recording of a hippocampal neuron while the rat runs back and forth in a
% linear maze.
%
% Data:
% expTime: time axis for entire experiment (in seconds at 1 ms resolution)
% ratPosition: rat's position (in cm) at each time point in expTime
% spikeTimes: time at which each recorded action potential occurred

%% load data
% make sure placeCellData.mat is in your path
load placeCellData.mat
nr = 4;
nc = 3;

%% plot the rat's position over time
main = figure('position',[50 50 1200 600]);
subplot(nr,1,1)
plot(expTime,ratPosition);
hold on;
xlabel('Time [s]');	
ylabel('Position [cm]')
title('Fig. 1: Rat position vs. time');

%% Plot spikes on top of position trace.
% Make a binary variable that is size 177761 x 1 that indicates when a
% spike occurred. Name that variable 'spikeTrain'. Use the space provided
% below. Hint: use  the variables 'expTime' and 'spikeTimes'

%ANSWER
spikeTrain = ismember(expTime,spikeTimes);
%

% Find the index of each spike and name that variable 'spikeIndex' (220 x 1).

%ANSWER
spikeIndex = find(spikeTrain);
%

% We then use that index to plot a dot of the rat's position at the time
% of that spike.
hp=plot(expTime(spikeIndex),ratPosition(spikeIndex),'r.');
set(hp,'MarkerSize',10);    % make dots bigger

%%QUESTION: When does the cell fire? Is it just a place cell?




%% Occupancy normalized histogram
% We want to visualize the probability of the cell firing as a function of
% position along the maze (ignoring for now the directionality issue).
% Because the rat is moving, it potentially spends more or less time in
% each spatial bin, so we need to normalize by the amount of time he spends
% in each bin.

positionBins=0:10:100;

% Using the positionBins indicated above, make a histogram of positions 
% where we got spikes. NOTE: You need to both plot a histogram and
% create a variable containing the spike-counts per bin so that you can create
% the normalized histogram below. The general way to do this is to use
% 'hist' to bin the data in a variable, then use 'bar' to create the plot:
% countsPerBin = hist(yData,positionBins);
% bar(positionBins,countsPerBin);

subplot(nr,nc,nc+1)

%ANSWER
spikeHist = hist(ratPosition(spikeIndex),positionBins);
bar(positionBins,spikeHist);
%

xlabel('Position [cm]')			%Label the axes.
ylabel('Spike count')
title('Spike histogram');
set(gca,'XTick',-50:50:150);  % Easier to see

% Using the positionBins indicated above, make a histogram of the occupancy
% times in seconds. Think carefully about what you are binning here. If,
% for example, a given position bin contains 100 counts (from the variable
% ratPosition), how man seconds did the rat spend in that position bin?
% As a reality check, you can see from fig. 1 that the entire experiment
% lasted just shy of 180 seconds. If you have calculated the occupancy
% times in seconds, then the sum of all occupancy bins should add up to the
% total length of the experiment.
subplot(nr,nc,nc+2)

%ANSWER
occupancyHist = hist(ratPosition,positionBins)/1000;
bar(positionBins,occupancyHist);
%

xlabel('Position [cm]')			%Label the axes.
ylabel('Time in spatial bin (s)')
title('Position histogram');

% Now make a histogram of the positions where spikes occurred that is 
% normalized by the occupancy time in each bin.
subplot(nr,nc,nc+3)

%ANSWER
bar(positionBins,spikeHist./occupancyHist);
%

xlabel('Position [cm]')			%Label the axes.
ylabel('Occupancy normalized counts (spikes/s)')
title('Occupancy normalized histogram');

% Compare the histogram in the lower left panel ('Spike histogram') with
% the one on the lower right ('Occupancy normalized histogram'). Are there
% any differences? 

%% Chapter 9, Model #1
% We want to fit a model that will predict the cell's spike counts in each
% bin as a function of its position along the track. The natural model 
% is the Poisson, where we express the mean rate as a function of time in
% terms of the covariates: lambda_t = beta0 + beta_1(position_t)

% Fit Poisson Model to the spike train data using the rat's position as a 
% predictor. Fill in the inputs below. See help on function 'glmfit'. 

%ANSWER
[b1,dev1,stats1] = glmfit(ratPosition,spikeTrain,'poisson','log');
%

%QUESTION
% What are each of these output terms (b1,dev1,stats1)?
%




%re-plot occupancy norm. hist.
subplot(nr,nc,2*nc+1)
bar(positionBins,spikeHist./occupancyHist);
hold on;
%Plot the model.
plot(positionBins,exp(b1(1)+b1(2)*positionBins)*1000,'r');
xlabel('Position [cm]')				%Label the axes.
ylabel('Occupancy normalized counts (spikes/s)')
title('Model 1: Position only covariate');

%% Improve model fit by adding a squared term for position: Model #2

%ANSWER
[b2,dev2,stats2] = glmfit([ratPosition, ratPosition.^2],spikeTrain,'poisson','log');
%

% look at the fit
subplot(nr,nc,2*nc+2)
bar(positionBins,spikeHist./occupancyHist);
hold on;                        	
plot(positionBins,exp(b2(1)+b2(2)*positionBins+b2(3)*positionBins.^2)*1000,'r');
xlabel('Position [cm]')				
ylabel('Occupancy normalized counts (spikes/s)')
title('Model 2: Position and Position-squared as covariates');
%Compute maximum likelihood estimates of:
mu=-b2(2)/2/b2(3);                  %...place field center,
sigma=sqrt(-1/(2*b2(3)));           %...place field size,
alpha=exp(b2(1)-b2(2)^2/4/b2(3));   %...max firing rate.

%% Let's look at how our model does on the raw data: raw residuals
cumResid = cumsum(stats2.resid);

subplot(nr,1,4);
yyaxis left
plot(expTime,cumResid);
xlabel('Time (s)');
ylabel('Cumulative residuals');

yyaxis right
plot(expTime,ratPosition);
ylabel('Position (cm)');

% Q1: Is there any relationship between the residuals of our model and the
% direction of motion of the rat

%ANSWER
% residuals get bigger as rat moves in positive direction and smaller
% as he moves in the negative direction.
%

%% Include direction of motion in the model: Model #3

% So we need a covariate for direction. Let's start with a simple
% indicator variable: 1 if rat is moving in positive direction, 0 otherwise

%ANSWER
ratDirection = [0; diff(ratPosition) > 0];
%

% Now we just throw this into the model as another covariate:
[b3,dev3,stats3] = glmfit([ratPosition,ratPosition.^2,ratDirection],spikeTrain,'poisson','log');

% Is the directional coefficient statistically significant? Check the p
% value in our stats output variable for each predictor.

%ANSWER
pBeta3 = stats3.p(4)<.05;    % darn tootin'
%

% and now re-do our cumulative residuals plot: much better
cumResid = cumsum(stats3.resid);
subplot(nr,1,4);
hold on
plot(expTime,cumResid,'k');
xlabel('Time (s)');
ylabel('Cumulative residuals');

%% Plot occupancy normalized histogram for each direction of motion separately

%EXTRA CREDIT QUESTION
% Why might it be useful to fit separate models for each direction of
% motion. Think about other predictors that we don't have access to. Are
% there any more predictors we could obtain from the data in our workspace?
% See if you can obtain any differential statistics of the animal's
% behavior in the forward and reverse directions. Examine how we can split
% model fit in forward and reverse directions below.






spikeTrainUp = spikeTrain & ratDirection;
spikeTrainDown = spikeTrain & ~ratDirection;
spikeIndexUp=find(spikeTrainUp);%Determine index of each spike.
spikeIndexDown=find(spikeTrainDown);%Determine index of each spike.
positionBins=0:10:100;

% Histogram of positions where we got spikes.
spikeHistUp=hist(ratPosition(spikeIndexUp),positionBins);
spikeHistDown=hist(ratPosition(spikeIndexDown),positionBins);

occupancyHist = hist(ratPosition,positionBins) .* (0.001/2);
subplot(nr,nc,[2*nc+3]);
hB = bar(positionBins,[(spikeHistUp./occupancyHist)',(spikeHistDown./occupancyHist)']);
hB(2).FaceColor = [1 0 0];
hold on;
xlabel('Position (cm)')			%Label the axes.
ylabel('Occupancy normalized counts (spikes/s)')
legend('Up','Down');
title('Occupancy normalized histograms for each direction');

%% Visualize the fit for each direction using glmval
positionBins=(0:10:100)';
nBins = size(positionBins);

% Evaluate our direction model in direction 0 (position decreases):
% glmval returns both the mean and the upper and lower CIs
[lambdaDown, CIloDown, CIhiDown] = glmval(b3,[positionBins,positionBins.^2,zeros(nBins)],...
    'log',stats3);

% Evaluate our direction model in direction 1 (position increases):
[lambdaUp, CIloUp, CIhiUp] = glmval(b3,[positionBins,positionBins.^2,ones(nBins)],...
    'log',stats3);

% plot the results
errorbar(positionBins,lambdaUp.*1000,CIloUp.*1000,CIhiUp.*1000,'b');
errorbar(positionBins,lambdaDown.*1000,CIloDown.*1000,CIhiDown.*1000,'r');
xlabel('Position (cm)');
ylabel('Firing rate (spikes/s)');


%%  Measures of goodness of fit

% "There is not a single procedure for measuring goodness-of-fit; instead
% there are many tools that, taken together, can provide a broad
% perspective on the strenghts and weaknesses of a set of models." p. 280

%% Method 1: Comparing Akaike's Information Criterions (AIC) values
% AIC is a form of "penalized likelihood" measure: we first compute
% -2*log(likelihood) of the data given the model (will ebe small for good
% models) and then add a penalty term "2*p," where p is the number of
% parameters in the model. 

%QUESTION
% Why do we penalize for the number of parameters in the model?





% Recall that what we are actually predicting is the Poisson rate
% parameter, lambda. For model #1 (only x as covariate)
lambda1 = exp(b1(1)+ b1(2)*ratPosition);          %Poisson rate function,
loglikelihood1 = sum(log(poisspdf(spikeTrain,lambda1)));     %log likelihood

%Calculate AIC for Model 1
AIC1 = -2 * loglikelihood1 + (2*2);                           %AIC for Model 1 (2 free params)




%Calculate AIC for Model 2
AIC2 = -2*sum(log(poisspdf(spikeTrain, ...
              exp(b2(1) + b2(2)*ratPosition + b2(3)*ratPosition.^2)))) + (2*3);
          
          
          

dAIC=AIC1-AIC2;		%Difference in AIC between Models 1 and 2; Your answer should be 636.0145.




% QUESTION
% What does this number mean? 


% ANSWER
% In general, we just pick the model with the
% lower AIC value, but we can also think of it in terms of the parameter
% penalty: the fact that model 2 has an AIC of ~636 less than the AIC of
% model 1 means that model 2 would still be preferable to model 1 even if
% it had 636/2 (= 318) more parameters than it actually does.



% NOTE: We can also more easily calculate AIC from the deviance 
% (The deviance is a generalization of the residual sum of squares for all 
% exponential family distributions. Sum of squares is only appropriate for 
% gaussian distributions.)

% deltaAIC = (Dev1 + 2*p1) - (Dev2 + 2*p2)
alt_dAIC = (dev1 + 2*2) - (dev2 + 2*3);     % compare with dAIC above

%% Method 3: Confidence intervals on model parameters

%Compute 95% CI for parameters of Model 1.
CI1 =[b1 - 2*stats1.se, b1 + 2*stats1.se];
eCI1 = exp(CI1);	%Exponentiate Model 1 CIs.

%Compute 95% CI for parameters of Model 2.
CI2 = [b2 - 2*stats2.se, b2 + 2*stats2.se];
eCI2 = exp(CI2);

% Wald test: hypothesis test for whether a parameter is significantly
% different from 0
pBeta2 = stats2.p(3);	%Significance level of Model 2 additional parameter.


%% Comparing model #3 (with direction term) vs. model #2

% Calculate the dAIC between Model 2 and Model 3

%ANSWER
dAIC = (dev2 + 2*3) - (dev3 + 2*4);
%

p2vs3 = 1 - chi2cdf(dev2-dev3,1);       % Wald test

%For model 3, compute 95% CI for last parameter and find the significance
%level using the Wald test.

%ANSWER
CI3 = [b3(4)-2*stats3.se(4), b3(4)+2*stats3.se(4)];

pBeta3 = stats3.p(4);	%... and significance level.
%

%QUESTION
%What do these results tell us about our three models?




