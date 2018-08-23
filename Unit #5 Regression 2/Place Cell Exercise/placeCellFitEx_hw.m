% placeCellFitEx_hw.m
% 
% This exercise has been adapted from:
% Chapter 9, of "Case studies in neural data analysis" by Kramer & Eden 2016
% 
% Homework instructions:
%
% Collaboration and consultation of outside sources. The goal of this
% exercise is to give you practice fitting a regression model in MATLAB.
% The instructions and supplied code below will lead you through this in a
% straightforward and gentle manner. Therefore, we ask that you do not
% collaborate with any classmates on this exercise prior to class
% discussion. Neither should you seek out answers from published materials
% (i.e. Kramer & Eden 2016) or online code repositories. This is a
% self-contained exercise, so you should just work through it on you own.
%
% What to do: Each section begins with a '%%'. Read through the comments
% and follow the instructions provided. In some cases you will be asked to
% answer a question, clearly indicated by 'QUESTION'--please do so using
% comments (i.e. by starting each line with a '%'). In other cases, you be
% asked to supply missing code, indicated by 'TODO'. Once you have supplied
% the required code, you can execute that section by mouse-clicking in that
% section (The block will turn yellow.) and then simultaneously hitting the
% 'ctrl' and 'enter' keys (PC) or 'command' and 'enter' keys (Mac).
%
% The 1st two sections of code don't require you to do anything, but you
% still need to execute them in order to load the data and display an
% initial plot.
%
% Original source of exercise:
% Mark A. Kramer and Uri T. Eden, "Case Studies in Neural Data Analysis",
% MIT Press, 2016, Chapter 9.
% https://mitpress.mit.edu/books/case-studies-neural-data-analysis
%
% Adapted by RTB, 9 Dec. 2016
% Developed for homework by LD and RTB, March-April 2017

%% Concepts covered:
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

% TODO: Make a binary variable that is size 177761 x 1 that indicates when a
% spike occurred. Name that variable 'spikeTrain'. Use the space provided
% below. Hint: use  the variables 'expTime' and 'spikeTimes'


spikeTrain = zeros(size(expTime));
tic; spikeTrain = hist(spikeTimes,expTime); toc

% TODO: Using spikeTrain, find the index of each spike and name that variable
% 'spikeIndex' (220 x 1).


spikeIndex = find(spikeTrain);


% We then use that index to plot a dot of the rat's position at the time
% of that spike.
hp=plot(expTime(logical(spikeTrain)),ratPosition(logical(spikeTrain)),'r.');
set(hp,'MarkerSize',10);    % make dots bigger

% QUESTION: When does the cell fire? Is it just a place cell?





%% Occupancy normalized histogram
% We want to visualize the probability of the cell firing as a function of
% position along the maze (ignoring for now the directionality issue).
% Because the rat is moving, it potentially spends more or less time in
% each spatial bin, so we need to normalize by the amount of time he spends
% in each bin.

positionBins=0:10:100;

% TODO: Using the positionBins indicated above, make a histogram of positions 
% where we got spikes. NOTE: You need to both plot a histogram and
% create a variable containing the spike-counts per bin so that you can create
% the normalized histogram below. The general way to do this is to use
% 'hist' to bin the data in a variable, then use 'bar' to create the plot.
% Name this variable 'spikeHist' (it will be used later).
% spikeHist = hist(yData,positionBins);
% bar(positionBins,spikeHist);
subplot(nr,nc,nc+1)
spikeHist = hist(ratPosition(spikeIndex),positionBins);
bar(positionBins,spikeHist);
xlabel('Position [cm]')			%Label the axes.
ylabel('Spike count')
title('Spike histogram');
set(gca,'XTick',-50:50:150);  % Easier to see

% TODO: Using the positionBins indicated above, make a histogram of the
% occupancy times in seconds. Think carefully about what you are binning
% here. If, for example, a given position bin contains 100 counts (from the
% variable ratPosition), how man seconds did the rat spend in that position
% bin? As a reality check, you can see from fig. 1 that the entire
% experiment lasted just shy of 180 seconds. If you have calculated the
% occupancy times in seconds, then the sum of all occupancy bins should add
% up to the total length of the experiment. Name this variable
% 'occupancyHist'
subplot(nr,nc,nc+2)
occupancyHist = hist(ratPosition,positionBins) ./ 1000;
bar(positionBins,occupancyHist);
xlabel('Position [cm]')			%Label the axes.
ylabel('Time in spatial bin (s)')
title('Position histogram');

% TODO: Now make a histogram of the positions where spikes occurred that is 
% adjusted for the rat's occupancy time in each position bin.
subplot(nr,nc,nc+3)
bar(positionBins,spikeHist ./ occupancyHist);
xlabel('Position [cm]')			%Label the axes.
ylabel('Occupancy normalized counts (spikes/s)')
title('Occupancy normalized histogram');

% QUESTION: Compare the histogram in the lower left panel ('Spike
% histogram') with the one on the lower right ('Occupancy normalized
% histogram'). Are there any differences?

%% Chapter 9, Model #1

% We want to fit a model that will predict the cell's spike counts in each
% bin as a function of its position along the track. The natural model 
% is the Poisson, where we express the mean rate as a function of time in
% terms of the covariates: lambda_t = beta0 + beta_1(position_t)
% However, as we discussed in lecture, the right side of our equation is
% not bounded and can assume negative values, whereas spike rates cannot go
% below 0. So the trick is to use a so-called 'link function' to transform
% our dependent variable--in this case we use the natural logarithm
% ('log'), so that we are actually fitting log(lambda_t). The generalized
% linear model makes this very easy. We don't even have to overtly take
% the log of our dependent variable; we just treat it like ordinary linear
% regression, but, in addition, specifiy the appropriate probability
% distribution (in our case, 'poisson') and the appropriate link function
% ('log'). But we must remember that we are still really fitting
% log(lambda_t) if we are to interpret our beta coefficients properly.

% TODO: Fit a Poisson Model to the spike train data using the rat's position as a 
% predictor. Fill in the inputs below. See help on function 'glmfit'. 


[b1,dev1,stats1] = glmfit(ratPosition,spikeTrain,'poisson','link','log');


%QUESTION: What are each of these output terms (b1,dev1,stats1)?

% Interpreting our beta coefficients is a bit trickier now, because they
% are predicting log(lambda_t), not lambda_t directly, and it is the latter
% in which we are interested. So in order to interpret the coefficients in
% a straightforward way, we need to 'undo' the natural logarithm.

% QUESTION: What is our predicted firing rate when the rat is at position
% 0? (HINT: Write down the model!)


%re-plot occupancy normalized histogram
subplot(nr,nc,2*nc+1)
bar(positionBins,spikeHist./occupancyHist);
hold on;
%Plot the model.
plot(positionBins,exp(b1(1)+b1(2)*positionBins)*1000,'r');
xlabel('Position [cm]')				%Label the axes.
ylabel('Occupancy normalized counts (spikes/s)')
title('Model 1: Position only covariate');

%% Model #2

% TODO: Improve model fit by ADDING a squared term for position.
% Hint: 'b2' should contain 3 elements

[b2,dev2,stats2] = glmfit([ratPosition,ratPosition.^2],spikeTrain,'poisson');


% Look at the fit
subplot(nr,nc,2*nc+2)
bar(positionBins,spikeHist./occupancyHist);
hold on;                        	
plot(positionBins,exp(b2(1)+b2(2)*positionBins+b2(3)*positionBins.^2)*1000,'r');
xlabel('Position [cm]')				
ylabel('Occupancy normalized counts (spikes/s)')
title('Model 2: Position and Position-squared as covariates');

% QUESTION: What kind of statistical distribution does our model resemble?



%% Re-cast the model for easier interpretation of beta coefficients

% We notice that model #2 looks sort of like a Gaussian (ans. to Q above!)
% And if we compare the model to the formula for a Gaussian, we notice that
% we can transform one into the other:
% Gaussian: lambda_t = alpha * exp((ratPosition - mu).^2 / (2*sigma.^2))
% Model 2: lambda_t = exp(beta0 + beta1*ratPosition + beta2*ratPosition^2);
% We've essentially replaced our 3 beta terms with 3 new terms (alpha, mu
% and sigma) that correspond to more intuitive concepts. With a little
% algebra, we can express the 3 Gaussian parameters in terms of our beta
% parameters.

% (As a cool parenthetical, it turns out from statistical theory that the
% maximum likelihood estimate (MLE) of any function of the fit parameters
% is the just the same function applied to the MLE parameters.)

%Compute maximum likelihood estimates of:
mu=-b2(2)/2/b2(3);                  %...place field center,
sigma=sqrt(-1/(2*b2(3)));           %...place field size,
alpha=exp(b2(1)-b2(2)^2/4/b2(3));   %...max firing rate.

%% Analysis of residuals

% Residuals tell us, on a point-by-point basis, the difference between the
% data and the predictions of our model. 'glmfit' gives us these values for
% free, and they are a valuable resource for evaluating deficiencies in our
% model. One type of residual that is particularly useful for spike data is
% the cumulative raw residual, which is just the sum of the residuals up to
% each point in time:
cumResid = cumsum(stats2.resid);

% Superimpose cumulative residuals on ratPosition over time
subplot(nr,1,4);
yyaxis left
plot(expTime,cumResid);
xlabel('Time (s)');
ylabel('Cumulative residuals');

yyaxis right
plot(expTime,ratPosition);
ylabel('Position (cm)');

% QUESTION: Is there any relationship between the residuals of our model
% and the direction of motion of the rat?

% QUESTION: What do you think the source of this relationship is?

% QUESTION: If we had the correct model, what should the cumulative
% residuals look like in a similar plot?

%% Include direction of motion in the model: Model #3

% TODO: To provide a covariate for direction, create an indicator variable,
% ratDirection, in which each bin contains a 1 if rat is moving in the
% positive direction, and 0 otherwise. (HINT: Your direction variable needs
% to be the same size as ratPosition and spikeTrain.)


ratDirection = [0;diff(ratPosition)>0];


% Now we just throw this into the model as another covariate:
[b3,dev3,stats3] = glmfit([ratPosition,ratPosition.^2,ratDirection],spikeTrain,'poisson','log');

% QUESTION: Is the directional coefficient statistically significant? Check
% the p-value in our stats output variable for each predictor. What is the
% relevant p-value for ratDirection?



% and now re-do our cumulative residuals plot: much better
cumResid = cumsum(stats3.resid);
subplot(nr,1,4);
hold on
plot(expTime,cumResid,'k');
xlabel('Time (s)');
ylabel('Cumulative residuals');

%% Plot occupancy normalized histogram for each direction of motion separately

% QUESTION (EXTRA CREDIT)
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
%subplot(nr,nc,[2*nc+3]);
figure;
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
% perspective on the strenghts and weaknesses of a set of models." 
% - Kramer & Eden 2016, p. 280

%% Method 1: Comparing Akaike's Information Criterion (AIC) values

% AIC is a form of "penalized likelihood" measure: we first compute
% -2*log(likelihood) of the data given the model (will be small for good
% models) and then add a penalty term "2*p," where p is the number of
% parameters in the model. 

%QUESTION
% Why do we penalize for the number of parameters in the model?





% Recall that what we are actually predicting is the Poisson rate
% parameter, lambda. For model #1 (only ratPosition as covariate)
% We first use our model to calculate the Poisson rate function
lambda1 = exp(b1(1)+ b1(2)*ratPosition);
% Then we ask, for each timebin, how likely is the data (spikeTrain) given
% our model's predicted rate parameter (lambda1). We take the log of the
% likelihood in each timebin and sum them all up.
loglikelihood1 = sum(log(poisspdf(spikeTrain,lambda1)));     %log likelihood

% TODO: Calculate AIC for Model 1
AIC1 = -2*loglikelihood1 + 2*2



% TODO: Calculate AIC for Model 2
lambda2 = exp(b2(1)+ b2(2)*ratPosition + b2(3)*ratPosition.^2);
loglikelihood2 = sum(log(poisspdf(spikeTrain,lambda2)));
AIC2 = -2*loglikelihood2 + 2*3 
%Difference in AIC between Models 1 and 2; Your answer should be 636 (and
%change)
dAIC=AIC1-AIC2;

% QUESTION: What does this number mean? 





% NOTE: We can also more easily calculate AIC from the deviance 
% (The deviance is a generalization of the residual sum of squares for all 
% exponential family distributions. Sum of squares is only appropriate for 
% Gaussian distributions.)
alt_dAIC = (dev1 + 2*2) - (dev2 + 2*3);     % compare with dAIC above

%% Method 2: Confidence intervals on model parameters

%Compute 95% CI for parameters of Model 1.
CI1 =[b1 - 2*stats1.se, b1 + 2*stats1.se];
eCI1 = exp(CI1);	%Exponentiate Model 1 CIs.

%Compute 95% CI for parameters of Model 2.
CI2 = [b2 - 2*stats2.se, b2 + 2*stats2.se];
eCI2 = exp(CI2);

% QUESTION: Why do we multiply the parameter's standard error by 2 to get a
% 95% CI?

% QUESTION: How can we use these confidence intervals to perform hypothesis
% testing on our parameters?


% MATLAB provides a p-value for each beta parameter. For example, the 
% significance level of Model 2's additional parameter is:
pBeta2 = stats2.p(3);

% QUESTION: What is the correpsonding null hypothesis?

% QUESTION: What is the relationship between the p-value we obtained
% for the additional parameter with an Nth% confidence interval?



%% Comparing model #3 (with direction term) vs. model #2

% TODO: Calculate the dAIC between Model 2 and Model 3





% TODO: For model 3, compute 95% CI for ratDirection parameter and find the
% significance level.





%QUESTION: What do these results tell us about our three models?


