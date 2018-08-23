% placeCellFitEx.m
% 
% Chapter 9, of "Case studies in neural data analysis" by Kramer & Eden 2016
%
% modified by RTB 9 Dec. 2016

% Concepts covered:
% 1. Working with spike data
% 2. Occupancy normalized histogram for place fields
% 3. Using glmfit to form a Poisson point-process model
% 4. Model selection through analysis of residuals
% 5. Model comparison through measures of goodness of fit: AIC,
%       Chi-square, parameter CIs, Kolmogorov-Smirnov plots
%
% Recording of a hippocampal neuron while a rat runs back and forth in a
% linear maze.
%
% Data variables:
% expTime: time axis for entire experiment (in seconds at 1 ms resolution)
% ratPosition: rat's position (in cm) at each time point in expTime
% spikeTimes: time at which each recorded action potential occurred

% In-class feedback on 17 April 2017:
% 1.	optional bonus question unclear
% 2.	more background on topics like AIC
% 3.	set office hours for those who need help


%% load data
%cd 'C:\usr\rick\doc\Committees\PIN\PIN Director\Courses\Stats\Eden Data\Chapter9'
load placeCellData.mat

%% plot the rat's position over time
figure, plot(expTime,ratPosition);
hold on;
xlabel('Time [s]');	
ylabel('Position [cm]')
title('Fig. 9.1: Rat position vs. time');
%% Plot spikes on top of position trace.
% This is a tad tricky. You actually want to plot the rat's POSITION at the
% time the spike occurred.

% First step is to create a vector the same size as expTime that has a 1 at
% each time-bin that contains a spike and 0's elsewhere. There are many
% ways to do this.
% A clever trick is to use hist by explicitly specifying bins
spikeTrain = hist(spikeTimes,expTime)';
% Now determine index of each spike . . .
spikeIndex=find(spikeTrain);
% . . . then use that index to plot a dot of the rat's position at the time
% of that spike.
hp=plot(expTime(spikeIndex),ratPosition(spikeIndex),'r.');
set(hp,'MarkerSize',10);    % make dots bigger

% When does the cell fire? Does it 'care' only about place?

%% Occupancy normalized histogram
% We want to visualize the probability of the cell firing as a function of
% position along the maze (ignoring for now the directionality issue).
% Because the rat is moving, it potentially spends more or less time in
% each spatial bin, so we need to normalize by the amount of time he spends
% in each bin (= occupancy).
% [NOTE: It would be good to have students plot both the raw spike
% histogram and compare it with the occupancy normalized version.]
positionBins=0:10:100;
% Histogram of positions where we got spikes.
spikeHist=hist(ratPosition(spikeIndex),positionBins);
% Histogram of the occupancy times in seconds
occupancyHist = hist(ratPosition,positionBins) .* 0.001;
figure;
bar(positionBins,spikeHist./occupancyHist);
xlabel('Position [cm]')			%Label the axes.
ylabel('Occupancy normalized counts (spikes/s)')
title('Fig. 9.4: Occupancy normalized histogram');
%% Chapter 9, Model #1
% We want to fit a model that will predict the cell's spike counts in each
% bin as a function of its position along the track. The natural model 
% is the Poisson, where we express the mean spike rate as a function of time in
% terms of the covariates: lambda_t = beta0 + beta_1(position_t)

% One intuitive way to think about this (following Oram's lecture on
% logistic regression) is that the right side of our eqn. must be able to
% go both pos. & neg., whereas the right, being a spike rate (and, really,
% a spike count), can only be a positive integer. So we take the log of
% this, which now allows us to include negative numbers.

%Fit Poisson Model to the spike train data.
[b1,dev1,stats1] = glmfit(ratPosition,spikeTrain,'poisson','log');
%re-plot occupancy norm. hist.
figure, bar(positionBins,spikeHist./occupancyHist);
hold on;
%Plot the model.
plot(positionBins,exp(b1(1)+b1(2)*positionBins)*1000,'r');
xlabel('Position [cm]')				%Label the axes.
ylabel('Occupancy normalized counts (spikes/s)')
title('Model 1: Position only covariate');

%% Improve model fit by adding a squared term for position: Model #2

[b2,dev2,stats2] = glmfit([ratPosition, ratPosition.^2],spikeTrain,'poisson','log');

% Look at the fit
figure, bar(positionBins,spikeHist./occupancyHist);
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
mu=-b2(2)/2/b2(3);                  %...place field center (63.16)
sigma=sqrt(-1/(2*b2(3)));           %...place field size, (9.57)
alpha=exp(b2(1)-b2(2)^2/4/b2(3));   %...max firing rate (0.0113)

%% Gaussian fit directly with unconstrained nonlinear optimization

OPTIONS = optimset('Display','off','TolX',0.001);
% Look at our occupancy normalized histogram to generate guesses
q0 = [20/1000,30,30];   % reasonable guesses for alpha, mu, sigma
%q0 = [0,0,0];           % BAD guesses
%qFit = fminsearch(@(q)fitFunGauss(q,ds),q0,OPTIONS);
qFit = fminsearch(@(q)fitFunGauss2(q,ratPosition,spikeTrain),q0,OPTIONS);

% Compare qFit with alpha, mu and sigma calculated above. Pretty good.
% But to see the down side of this approach, try making initial guesses
% that are less well guided by the histogram. Do we always converge to the
% correct answer?

%% Analysis of residuals

% Residuals tell us, on a point-by-point basis, the difference between the
% data and the predictions of our model. 'glmfit' gives us these values for
% free, so we may as well make use of them. They are a valuable resource
% for evaluating deficiencies in our model. One type of residual that is
% particularly valuable for spike data is the cumulative raw residual,
% which is just the sum of the residuals up to each point in time:
cumResid = cumsum(stats2.resid);

figure, subplot(2,1,1);
plot(expTime,cumResid);
ax = axis;  % for later comparison with directional model
xlabel('Time (s)');
ylabel('Cumulative residuals');
subplot(2,1,2);
plot(expTime,ratPosition);
xlabel('Time (s)');
ylabel('Position (cm)');

%% Superimpose residuals on position trace
% There seems to be a relationship between the cumulative residuals and the
% rat's position along the track. To better see this, let's superimpose
% them directly using a double-y axis plot. You may need to increase the
% size of this figure to more clearly see the relationships.

figure
yyaxis left
plot(expTime,cumResid);
xlabel('Time (s)');
ylabel('Cumulative residuals');

yyaxis right
plot(expTime,ratPosition);
ylabel('Position (cm)');

hold on
hp=plot(expTime(spikeIndex),ratPosition(spikeIndex),'k.');
set(hp,'MarkerSize',10);    % make dots bigger

% QUESTION: What is the pattern here? When do the cumulative residuals tend to go
% up and when do they tend to go down?


% QUESTION: If we had the correct model, what should the cumulative
% residuals look like in a similar plot?



%% Include direction of motion in the model: Model #3
% 

% So we need a covariate for direction. Let's start with a simple
% indicator variable: 1 if rat is moving in positive direction, 0 otherwise
ratDirection = [0; diff(ratPosition) > 0];

% Now we just throw this into the model as another covariate:
[b3,dev3,stats3] = glmfit([ratPosition,ratPosition.^2,ratDirection],spikeTrain,'poisson','log');

% Is the directional coefficient statistically significant?
pBeta3 = stats3.p(4);    % darn tootin'

% and now re-do our cumulative residuals plot: much better
cumResid = cumsum(stats3.resid);
figure, subplot(2,1,1);
plot(expTime,cumResid);
axis(ax);   % same axis as other residual plot, for comparison
xlabel('Time (s)');
ylabel('Cumulative residuals');
subplot(2,1,2);
plot(expTime,ratPosition);
xlabel('Time (s)');
ylabel('Position (cm)');

%% Plot occupancy normalized histogram for each direction of motion
spikeTrainUp = spikeTrain & ratDirection;
spikeTrainDown = spikeTrain & ~ratDirection;
spikeIndexUp=find(spikeTrainUp);%Determine index of each spike.
spikeIndexDown=find(spikeTrainDown);%Determine index of each spike.
positionBins=0:10:100;
% Histogram of positions where we got spikes.
spikeHistUp=hist(ratPosition(spikeIndexUp),positionBins);
spikeHistDown=hist(ratPosition(spikeIndexDown),positionBins);
% Histogram of the occupancy times in seconds. For the time being, we are
% assuming the rat moved at about the same speed in each direction. Does
% this look like a reasonable assumption? How could we make it more
% precise?
% Need to divide by 2 because we're only looking at one direction. This is
% a good teachable moment to make the students think about the nature of
% the occupancy normalized histogram, because the GLM will give the real
% values (i.e. they will look 2x too big if you don't fix the occupancy
% normalization procedure by dividing by 2).
occupancyHist = hist(ratPosition,positionBins) .* (0.001/2);
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

%return
%%  Measures of goodness of fit (work in progress)

% "There is not a single procedure for measuring goodness-of-fit; instead
% there are many tools that, taken together, can provide a broad
% perspective on the strenghts and weaknesses of a set of models." p. 280

%% Method 1: Comparing Akaike's Information Criterions (AIC) values
% AIC is a form of "penalized likelihood" measure: we first compute
% -2*log(likelihood) of the data given the model (will ebe small for good
% models) and then add a penalty term "2*p," where p is the number of
% parameters in the model.

% Recall that what we are actually predicting is the Poisson rate
% parameter, lambda. For model #1 (only x as covariate)
lambda1 = exp(b1(1)+ b1(2)*ratPosition);        %Poisson rate function,
LL1 = sum(log(poisspdf(spikeTrain,lambda1)));   %log likelihood
AIC1 = -2*LL1 + (2*length(b1));                 %AIC for Model 1 (2 free params)

%Compute the AIC for Model 2: position & position^2
AIC2 = -2*sum(log(poisspdf(spikeTrain, ...
              exp(b2(1) + b2(2)*ratPosition + b2(3)*ratPosition.^2)))) + (2*length(b2));

dAIC=AIC1-AIC2;		%Difference in AIC between Models 1 and 2; 636.0145

% What does this number mean? In general, we just pick the model with the
% lower AIC value, but we can also think of it in terms of the parameter
% penalty: the fact that model 2 has an AIC of ~636 less than the AIC of
% model 1 means that model 2 would still be preferable to model 1 even if
% it had 636/2 (= 318) more parameters than it actually does.

% NOTE: We can also more easily calculate AIC from the deviance 
% (The deviance is a generalization of the residual sum of squares.)
% deltaAIC = (Dev1 + 2*p1) - (Dev2 + 2*p2)
alt_dAIC = (dev1 + 4) - (dev2 + 6);     % compare with dAIC above

%% Method 2a: Chi-square test for nested models
% The AIC method does not tell us whether adding a parameter
% *significantly* improves our fit to the data. To do this, we can use a
% class of hypothesis test called "maximum likelihood ratio tests".

% For this test, H0 is that the additional parameters (i.e. those in the
% full model that are not contained in the reduced model) all = 0. The test
% statistic for the MLRT is equivalent to the difference in deviances:
% S = dev1 - dev2
% Under H0, this statistic follows a chi-square distribution with n2 - n1
% (i.e. diff. in # of parameters) degrees of freedom:
pValChiSq = 1-chi2cdf(dev1 - dev2,1);	%Compare Models 1 and 2, nested GLMs.

%% Method 2b: Sequential F-test for nested models:
% In general, the F statistic is a ratio of variances. In this case, the
% numerator is the additional variance accounted for by the full model (per
% additional free parameter), and the denominator is the variance not
% accounted for by the full model.
% That is: 
% Fnumerator = [(regression SS for higher degree model) – (regression SS for lower degree model)] /
% (difference in number of free parameters in the two models) = (regSS_full
% – regSS_red) / (dfF – dfR); But note that regSS = TotSS – resSS, so that
% the above can be written: = (resSS_red – resSS_full) / (dfF – dfR); 
% or, in Matlab, using the above variables: = (res1 – res2) / (length(qFit2) – length(qFit1));
% Fdenominator = residual MS for higher degree model = res2 / (n – length(qFit2))

%Fval = ((dev1-dev2) ./ (length(b2)-length(b1))) ./ (dev2 ./ (stats2.dfe));
resSSfull = sum(stats2.resid.^2);   % residual sum of squares for full model
resSSred = sum(stats1.resid.^2);    % residual sum of squares for reduced model
Fval = ((resSSred - resSSfull) ./ (length(b2)-length(b1))) ./ (resSSfull ./ (stats2.dfe));
pValFtest = 1 - fcdf(Fval, length(b2)-length(b1), stats2.dfe);

%% Method 3: Confidence intervals on model parameters

%Compute 95% CI for parameters of Model 1.
CI1 =[b1 - 2*stats1.se, b1 + 2*stats1.se];
eCI1 = exp(CI1);	%Exponentiate Model 1 CIs.

% Confidence intervals:
% b1(1) = [0.0004 to 0.0008], ~ 0.4 to 0.8 spikes per second at X = 0
% b1(2) = [1.0090 to 1.0171], ~ each 1 cm. increase in position leads to an
% increase in firing rate of between 0.9% and 1.7%.

%Compute 95% CI for parameters of Model 2.
CI2 = [b2 - 2*stats2.se, b2 + 2*stats2.se];
eCI2 = exp(CI2);

% Wald test: hypothesis test for whether a parameter is significantly
% different from 0
pBeta2 = stats2.p(3);	%Significance level of Model 2 additional parameter.

%% Method 4: Kolmogorov-Smirnov test
nSpikes = length(spikeTimes);				%Define # of spikes
%Evaluate Model 2.
lambda2 = exp(b2(1) + b2(2)*ratPosition + b2(3)*ratPosition.^2);

% Re-scale waiting times (see p. 286-7 for explanation)
Z = zeros(nSpikes,1);
Z(1) = sum(lambda2(1:spikeIndex(1)));	%1st rescaled waiting time.
for i = 2:nSpikes                             %... and the rest.
  Z(i) = sum(lambda2(spikeIndex(i-1):spikeIndex(i)));
end

%Compute empirical CDF from rescaled waiting times.
[eCDF, zVals] = ecdf(Z);

mCDF = 1-exp(-zVals);                       %Model CDF at z values.
figure, plot(mCDF,eCDF)						%Create KS-plot.
hold on
plot([0 1], [0 1]+1.36/sqrt(nSpikes),'k--')	%Upper confidence bound.
plot([0 1], [0 1]-1.36/sqrt(nSpikes),'k--')	%Lower confidence bound.
xlabel('Model CDF')                         %Label the axes.
ylabel('Empirical CDF')
axis([0,1,0,1]);
title('Fig. 9.8: KS plot of rescaled data for model 2');

%% Comparing model #3 (with direction term) vs. model #2
dAIC = (dev2 + 2*length(b2)) - (dev3 + 2*length(b3));
p2vs3 = 1 - chi2cdf(dev2-dev3,length(b3) - length(b2));       % Wald test

% Sequential F-test
resSSfull = sum(stats3.resid.^2);   % residual sum of squares for full model
resSSred = sum(stats2.resid.^2);    % residual sum of squares for reduced model
Fval = ((resSSred - resSSfull) ./ (length(b3)-length(b2))) ./ (resSSfull ./ (stats3.dfe));
pValFtest = 1 - fcdf(Fval, length(b3)-length(b2), stats3.dfe);

%For model 3, compute 95% CI for last parameter,
CIbeta3 = [b3(4)-2*stats3.se(4), b3(4)+2*stats3.se(4)];
p_beta3 = stats3.p(4);	%... and significance level.

% KS plot for model #3
lambda3 = exp(b3(1) + b3(2)*ratPosition + b3(3)*ratPosition.^2 + b3(4)*ratDirection);

% Re-scale the waiting times
Z = zeros(nSpikes,1);
Z(1) = sum(lambda3(1:spikeIndex(1)));	%1st rescaled waiting time.
for i=2:nSpikes							%... and the rest.
  Z(i)=sum(lambda3(spikeIndex(i-1):spikeIndex(i)));
end

[eCDF, zVals] = ecdf(Z);                    %Define empirical CDF,
mCDF = 1-exp(-zVals);                       %...and model CDF,
figure, plot(mCDF,eCDF)                     %...to create KS-plot.
hold on
plot([0 1], [0 1]+1.36/sqrt(nSpikes),'k')	%Upper confidence bound.
plot([0 1], [0 1]-1.36/sqrt(nSpikes),'k')	%Lower confidence bound.
xlabel('Model CDF'); ylabel('Empirical CDF');
axis([0,1,0,1]);
title('Figure 9.10: KS plot of rescaled data for model #3');

% Overall conclusion: Model #3 = pretty darned good.