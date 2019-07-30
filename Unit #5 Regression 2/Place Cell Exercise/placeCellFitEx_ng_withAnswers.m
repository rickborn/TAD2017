% placeCellFitEx_ng_withAnswers.m
% 
% This exercise has been adapted from:
% Chapter 9, of "Case studies in neural data analysis" by Kramer & Eden 2016
% 
% What to do: Login to Learning Catalytics (LC) and join the session for
% the module entitled "placeCellFit". You will answer a series of questions
% based on the guided programming below. Each section begins with a '%%'.
% Read through the comments and follow the instructions provided. In some
% cases you will be asked to answer a question, clearly indicated by
% 'QUESTION' and a corresponding 'Q#' that directs you to answer the
% relevant question in LC. In other cases, you be asked to supply missing
% code, indicated by 'TODO'. Once you have supplied the required code, you
% can execute that section by mouse-clicking in that section (The block
% will turn yellow.) and then simultaneously hitting the 'ctrl' and 'enter'
% keys (PC) or 'command' and 'enter' keys (Mac).
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
% Refined by RTB, September 2017

% % Concepts covered:
% 1. Working with spike data: times to indices
% 2. Occupancy normalized histogram for place fields
% 3. Using glmfit to form a Poisson point-process model
% 4. Model selection through analysis of residuals
% 5. Model comparison through measures of goodness of fit: AIC,
%       Chi-square, parameter CIs, Kolmogorov-Smirnov 
%
% Scenario: A rat runs back and forth in a linear maze, incentivized by a
% morsel of chocolate each time she reaches the end of the maze. While she
% is running we are recording the action potentials from a single neuron in
% the rat's hippocampus.
%
% Data: 
% - expTime: time axis for entire experiment
%   (in seconds at 1 ms resolution) 
% - ratPosition: rat's position (in cm) at each time point in expTime 
% - spikeTimes: time at which each recorded action potential occurred
%   (in seconds at 1 ms resolution)
%
% Scientific question: What aspects of the rat's behavior influence the
% neuron's tendency to fire an action potential?

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

% Look at the 3 experimental variables and make sure you understand what
% each one means and how the data are formatted.

% QUESTION (Q1): 
% What was the duration of the entire experiment, in seconds (rounded to
% nearest second)?

% QUESTION (Q2):
% How many action potentials did the hippocampal neuron fire during the
% entire experiment?

%% Plot spikes on top of position trace.

% Make a binary variable that is size 177761 x 1 that indicates when a
% spike occurred. Name that variable 'spikeTrain'. Use the space provided
% below. Hint: use  the variables 'expTime' and 'spikeTimes'

%ANSWER
spikeTrain = ismember(expTime,spikeTimes);
%
% Alternate ANSWER:
% spikeTrain = hist(spikeTimes,expTime);

% Find the index of each spike and name that variable 'spikeIndex' (220 x 1).
%ANSWER
spikeIndex = find(spikeTrain);

% We then use that index to plot a dot of the rat's position at the time
% of that spike.
hp=plot(expTime(spikeIndex),ratPosition(spikeIndex),'r.');
set(hp,'MarkerSize',10);    % make dots bigger

% QUESTION: When does the cell fire? Is it just a place cell?

%% Occupancy normalized histogram

% We want to visualize the probability of the cell firing as a function of
% position along the maze (ignoring for now the directionality issue).
% Because the rat may not be moving at a perfectly constant speed, it
% potentially spends more or less time in each spatial bin, so we need to
% normalize by the amount of time she spends in each bin.
positionBins = round(min(ratPosition)):10:round(max(ratPosition));

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

xlabel('Position [cm]')			%Label the axes.
ylabel('Spike count')
title('Spike histogram');
set(gca,'XTick',0:50:100);  % Easier to see
ax = axis;
axis([-10,110,ax(3),ax(4)]);

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

xlabel('Position [cm]')			%Label the axes.
ylabel('Time/bin (s)')
title('Position histogram');
set(gca,'XTick',0:50:100);  % Easier to see
ax = axis;
axis([-10,110,ax(3),ax(4)]);

% Now make a histogram of the positions where spikes occurred that is 
% normalized by the occupancy time in each bin.
subplot(nr,nc,nc+3)

%ANSWER
bar(positionBins,spikeHist./occupancyHist);

xlabel('Position [cm]')			%Label the axes.
ylabel('Occ. nl. rate (sp/s)')
title('Occupancy normalized histogram');
set(gca,'XTick',0:50:100);  % Easier to see
ax = axis;
axis([-10,110,ax(3),ax(4)]);

% Compare the histogram in the lower left panel ('Spike histogram') with
% the one on the lower right ('Occupancy normalized histogram'). Are there
% any differences? 

%% Model #1: Position as a covariate

% We want to fit a model that will predict the cell's spike counts in each
% bin as a function of its position along the track. The natural model 
% is the Poisson, where we express the mean rate as a function of time in
% terms of the covariates: lambda_t = beta0 + beta_1(position_t)

% Fit Poisson Model to the spike train data using the rat's position as a 
% predictor. Fill in the inputs below. See help on function 'glmfit'. 

%ANSWER
[b1,dev1,stats1] = glmfit(ratPosition,spikeTrain,'poisson','link','log');

%QUESTION
% What are each of these output terms (b1,dev1,stats1)?

% QUESTION (Q5): What is our predicted firing rate when the rat is at position
% 0? (HINT: Write down the model!)

% ANSWER: 
% log(E(Y|X)) = beta(1) + beta(2)'*X
% E(Y|X) = exp(intercept + beta(2)'*X)

% log(E(spikes | ratPosition)) = b1(1) + b1(2) .* ratPosition

% When the position is 0:
% E(Spikes | x = 0) = exp(intercept + b1(2)*0) = exp(intercept).
rate0 = exp(b1(1)) * 1000;  % 0.5879 spikes/sec

%re-plot occupancy norm. hist.
subplot(nr,nc,2*nc+1)
bar(positionBins,spikeHist./occupancyHist);
hold on;
%Plot the model.
plot(positionBins,exp(b1(1)+b1(2)*positionBins)*1000,'r');
xlabel('Position [cm]')				%Label the axes.
ylabel('Occ. nl. rate (sp/s)')
title('Model 1: Position only covariate');
set(gca,'XTick',0:50:100);  % Easier to see
ax = axis;
axis([-10,110,ax(3),ax(4)]);

%% Improve model fit by adding a squared term for position: Model #2

%ANSWER
[b2,dev2,stats2] = glmfit([ratPosition, ratPosition.^2],spikeTrain,'poisson','link','log');

% look at the fit
subplot(nr,nc,2*nc+2)
bar(positionBins,spikeHist./occupancyHist);
hold on;                        	
plot(positionBins,exp(b2(1)+b2(2)*positionBins+b2(3)*positionBins.^2)*1000,'r');
xlabel('Position [cm]')				
ylabel('Occ. nl. rate (sp/s)')
title('Model 2: Add Position^2 as covariate');
set(gca,'XTick',0:50:100);  % Easier to see
ax = axis;
axis([-10,110,ax(3),ax(4)]);

% QUESTION: What kind of statistical distribution does our model resemble?
% ANSWER: a Gaussian or normal distribution

%% Re-cast the model for easier interpretation of beta coefficients

% We notice that model #2 looks sort of like a Gaussian (ans. to Q6 above!)
% And if we compare the model to the formula for a Gaussian, we notice that
% we can transform one into the other:
%
% Gaussian: lambda_t = alpha * exp((ratPosition - mu).^2 / (2*sigma.^2))
% Model 2: lambda_t = exp(beta0 + beta1*ratPosition + beta2*ratPosition^2);
%
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

%% BONUS: Gaussian fit directly with unconstrained nonlinear optimization

% At this point, you might be asking yourself, "Why go through all of that
% algebraic gymnastics? Why not just fit a Gaussian model directly?" Well,
% we can do this using 'fminsearch'. What we need to do is
% write an objective function, call it 'fitFunGauss2', that will take in 3
% parameters and return some measure of error in the fit to the actual
% data. The 3 input parameters should be:
%
% 'q' : a vector of our 3 parameters (q(1)=alpha, q(2)=mu, q(3)=sigma)
% 'ratPosition' : our independent variable (x)
% 'spikeTrain' : our dependent variable (y)
%
% This function should first use the parameters to calculate values of
% lambda from ratPosition as predicted by the current model (i.e. the
% parameters in 'q'). Then it should compute the likelihood of the actual
% values ('spikeTrain') given the model's current estimates of the lambdas,
% and finally return the fit error as -log of the sum of these likelihoods:
%
%   -sum(log(poisspdf(y,lambda)))
%
% We then use 'fminsearch' to find the parameters, q, that minimize the
% value returned by 'fitFunGauss2'.

% THOUGHT QUESTION: Why do we have our function return *minus* the sum of
% the log likelihood? HINT: Think about what 'fminsearch' is doing and the
% overall goal, which is to maximize the likelihood of our data given the
% model.
%
% ANSWER: We are calculating a maximum likelihood estimate for our
% paramters by minimizing -log(likelihood) of the model given the data.

OPTIONS = optimset('Display','off','TolX',0.001);
% Look at our occupancy normalized histogram to generate guesses:
q0 = [20/1000,30,30];   % reasonable guesses for alpha, mu, sigma
qFit = fminsearch(@(q)fitFunGauss2(q,ratPosition,spikeTrain),q0,OPTIONS);

% Compare qFit with alpha, mu and sigma calculated above:
[qFit', [alpha;mu;sigma]]

% But to see the down side of this approach, try making initial guesses
% that are less well guided by the histogram (e.g. q0 = [0,0,0]). 
% Do we always converge to the correct answer? Ans. No

% THOUGHT QUESTION: Can you think of other benefits of using the GLM
% approach?
%
% ANSWER: glmfit is faster, we are guaranteed to get the right answer, and
% we automatically get all sorts of useful information, such as standard
% errors and confidence intervals.

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

% QUESTION (Q6): Describe the relationship between the residuals of our
% model and the position of the rat over time.

% ANSWER (A6): The residuals increase when the rat moves in the
% positive direction, don't change when the rat is not moving at either
% ends of the track and then decrease when the rat moves in the negative
% direction.

% QUESTION (Q7): Describe the source of this relationship. Think about what
% the current model "knows" (i.e. what predictor variables it contains) and
% what the cell is actually doing.

% ANSWER (A7): The firing rate of the cell depends on both position and
% direction of motion, but only position is accounted for in our current
% model. At the same positions, the cell fires a lot when the rat moves in
% one direction, but not at all when it moves in the opposite, so the model
% will fit the average. Thus the model's prediction is too low when the rat
% moves in the positive direction (which produces positive residuals) and
% too high when it moves in the negative direction (which produces negative
% residuals).

% QUESTION (Q8): If we had the "correct" model, what should the cumulative
% residuals look like in a similar plot?

% ANSWER (A8): They would be randomly scattered and close to 0. Note that
% even a "perfect" model would not produce residuals that were all zero,
% because you always have some irreducible noise that is, by definition,
% not fittable. What are the sources of such noise?

%% %% Model #3: Add direction of motion to the model

% So we need a covariate for direction. Let's start with a simple
% indicator variable: 1 if rat is moving in positive direction, 0 otherwise
ratDirection = [0; diff(ratPosition) > 0];

% Now we just throw this into the model as another covariate:
[b3,dev3,stats3] = glmfit([ratPosition,ratPosition.^2,ratDirection],spikeTrain,'poisson','link','log');

% or with 'fitglm' to get more diagnostics:
mdl3 = fitglm([ratPosition,ratPosition.^2,ratDirection],spikeTrain,'linear','Distribution','Poisson');

% QUESTION (Q9): Is the directional coefficient statistically significant?
% Check the p value in our stats output variable for each predictor.
% ANSWER (A9):
pBeta3 = stats3.p(4);    % darn tootin'

% and now re-do our cumulative residuals plot: much better
cumResid = cumsum(stats3.resid);
subplot(nr,1,4);
hold on
yyaxis left
plot(expTime,cumResid,'k');
xlabel('Time (s)');
ylabel('Cum. residuals');

% Add a line at 0 to help evaluate our residuals
ax = axis;
hl = line([ax(1),ax(2)],[0,0]);
set(hl,'Color','g','LineStyle','--');

%% Model parameter interpretation:

% Note that the log-link function changes how we interpret our beta
% coefficients. We already saw how to interpret beta0 above for Model #1.
% In the improved model (#3):
%
% Beta0:
% exp(b3(1))*1000 tells us that our predicted spike rate at position 0 is
% about 3e-10, or essentially 0. That's good.
%
% For the other parameters, we have to remember that:
%   exp(a+b) = exp(a) * exp(b)
%
% which means that a unit increase in the covariate has a multiplicative
% effect on the rate function:
%
% lambda_t = exp{b3(1) + b3(2)*position + b3(3)*position^2 + b3(4)*dir}
%
% so, for example, if we just consider the effect of position:
% exp(b3(2)) = 1.9915
% meaning that a 1 cm increase in the rat's position nearly doubles the
% cell's firing rate, whereas the effect of the position-squared term:
% exp(b3(3)) = 0.9946
% means we get about a 1% decrease for a unit change here. However, as the
% position increases, obviously position^2 increases much faster, so that
% gradually the 1% decrease overtakes the 200% increase and drives our
% predicted firing rate back down. Where in position this transition
% happens will determine the center of our place field, so it makes perfect
% sense that when we transformed our beta's into Gaussian parmeters (i.e.
% alpha, mu and sigma), the 'mu' parameter was determined only by the beta
% coefficients corresponding to the position and position^2 terms.
% (i.e. for Model #2, mu=-b2(2)/2/b2(3);)
%
% And finally, direction has a massive multiplier effect of 26.45, as it
% must in order to account for the fact that, at the exact same position,
% the cell fires a lot in one direction and not at all in the opposite.

% But NOTE: We can compute the same Gaussian description of the place field
% using our new parameters. The new parameter, direction, just acts as a
% multiplier of the same 3 transformed parameters:
mu = -b3(2) / (2*b3(3));                            %...place field center,
sigma = sqrt(-1/(2*b3(3)));                         %...place field size,
alpha = exp( b3(1) - ((b3(2)^2) / (4*b3(3))));      %...max firing rate.

%% Plot occupancy normalized histogram for each direction of motion separately

% EXTRA CREDIT QUESTION
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
ylabel('Occ. nl. rate (sp/s)')
legend('Up','Down','Location','Northwest');
title('Occ. nl. h-grams for each direction');
set(gca,'XTick',0:50:100);  % Easier to see
ax = axis;
axis([-10,110,ax(3),ax(4)]);

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
ylabel('FR (spikes/s)');

% QUESTION (Q10): Save this figure (fig. #1) as a jpeg and upload it to
% Learning Catalytics.

%%  Measures of goodness of fit

% "There is not a single procedure for measuring goodness-of-fit; instead
% there are many tools that, taken together, can provide a broad
% perspective on the strenghts and weaknesses of a set of models." 
% Kramer & Eden, p. 280

%% Method 1: Comparing Akaike's Information Criterions (AIC) values

% AIC is a form of "penalized likelihood" measure: we first compute
% -2*log(likelihood) of the data given the model (will be smaller for
% better models) and then add a penalty term "2*p," where p is the number
% of parameters in the model.

% QUESTION (Q11): Why do we include a penalty for the number of parameters
% in the model?

% ANSWER (A11): We can always improve our fit by adding more parameters to
% a model, so we need some sort of penalty that discourages such
% overfitting. AIC is a principled way of doing this. It uses information
% theory to balance between goodness-of-fit and model complexity. For more
% details, see:
%
% https://en.wikipedia.org/wiki/Akaike_information_criterion

% Recall that what we are actually predicting is the Poisson rate
% parameter, lambda. For model #1 (only position as covariate)
% This is model #1's prediction of the Poisson rate function:
lambda1 = exp(b1(1)+ b1(2)*ratPosition);

% Use 'poisspdf' to calculate the likelihood of our spiketrain given the
% model, then take the log of the probabilities and add them up:
logLikelihood1 = sum(log(poisspdf(spikeTrain,lambda1)));

% Calculate AIC for Model 1 (2 parameters)
AIC1 = -2 * logLikelihood1 + (2*2); 

% TODO: Calculate AIC for Model 2 (3 parameters)
lambda2 = exp(b2(1) + b2(2)*ratPosition + b2(3)*ratPosition.^2);
logLikelihood2 = sum(log(poisspdf(spikeTrain,lambda2)));
AIC2 = -2 * logLikelihood2 + (2*3);

% TODO: Calculate the difference in AIC values for Models 1 and 2:
dAIC = AIC1 - AIC2;

% QUESTION (Q12): What is the difference in AIC values between Model 1 and
% Model 2?
%
% ANSWER (A12): 636.

% QUESTION (Q13): What does this difference in AIC values mean? 
%
% ANSWER (A13): In general, we just pick the model with the lower AIC
% value, but we can also think of it in terms of the parameter penalty: the
% fact that model 2 has an AIC of ~636 less than the AIC of model 1 means
% that model 2 would still be preferable to model 1 even if it had 636/2 (=
% 318) more parameters than it actually does.

% NOTE: We can more easily calculate AIC from the deviance 
% (The deviance is a generalization of the residual sum of squares for all 
% exponential family distributions. Sum of squares is only appropriate for 
% gaussian distributions.)

% deltaAIC = (Dev1 + 2*p1) - (Dev2 + 2*p2)
alt_dAIC = (dev1 + 2*2) - (dev2 + 2*3);     % compare with dAIC above

%% Method 2: Confidence intervals on model parameters

% Recall that the linear regression framework gives us standard errors on
% our model parameters "for free." We can use these to compute confidence
% intervals, because maximum likelihood estimators are approximately normal
% for sufficiently large n.

% QUESTION (Q14): What is the standard error for beta0 (i.e. the
% y-intercept) in Model #1?
%
% ANSWER (A14): 0.15

% TODO: Calculate the 95% CI for the parameters of Model 1.
CI1 =[b1 - 1.96*stats1.se, b1 + 1.96*stats1.se];
% HINT: Recall that the linear model is a prediction of log(lambda). What
% do we do to interpret our beta's in terms of spike rate? 
eCI1 = exp(CI1) .* 1000;	%Exponentiate Model 1 CIs.

% QUESTION (Q15): For Model #1, what is the 95% confidence interval for the
% neuron's spiking rate at position x=0, in spikes per second?
%
% ANSWER (A15): [0.4401    0.7855]
%
% Conceptually, this is one of the more important questions in this
% exercise. Given the power of today's software, any moron can fit fancy
% regression models. The tricky part is in understanding the physical
% interpretation of our fit parameters. The first critical insight is that,
% for the Poisson point process model, we are predicting the logarithm of
% the rate function, not the rate function itself. So we need to
% exponentiate our betas to get to predictions about spike rate. The second
% thing you need to realize is that the predicted lambdas are for each time
% bin, and these are 1 ms wide. So the exponentiated CI of [0.0004 to
% 0.0008] is in spikes per millisecond and therefore need to be multipled
% by 1000 to get to spikes/sec.

% TODO: Compute the 95% CI for parameters of Model 2.
CI2 = [b2 - 1.96*stats2.se, b2 + 1.96*stats2.se];

% QUESTION (Q16): Based on your 95% CI, can we say that the
% position-squared term significantly (at alpha < 0.05) improves the model
% fit?
%
% ANSWER (A16): Yes. Since of 95% CI for the position-squared term [-0.0063
% -0.0046] does not include the null value of 0, we can reject H0.
%
% NOTE: The null value of 0 only applies to the un-exponentiated CI for our
% betas. Once we exponentiate, nothing can contain 0. If we are looking at
% the exponentiated values, then we have to use a null value of 1.

% We can also perform a direct signifcance test based on each parameter's
% maximum likelihood estimate and its standard error. This is called the
% 'Wald test', and we also get it for free with glmfit
%
% QUESTION (Q17): What is the p-value for the position-squared term in
% model #2?
%
% ANSWER (A17): 
pBeta2 = stats2.p(3);


%% Comparing model #3 (with direction term) vs. model #2

% QUESTION (Q18): What is the difference in AIC values between Model 2 and
% Model 3?
%
% ANSWER (A18): 
dAIC = (dev2 + 2*3) - (dev3 + 2*4);

% Model comparison with the Wald test: difference in deviances; the degrees
% of freedom is equal to the difference in #s of parameters in the 2
% models.
p2vs3 = 1 - chi2cdf(dev2-dev3,1);

% TODO: For model 3, compute 95% CI for direction parameter and find the
% significance level using the Wald test.
% ANSWER:
CI3 = [b3(4)-2*stats3.se(4), b3(4)+2*stats3.se(4)];
pBeta3 = stats3.p(4);	%... and significance level.

% QUESTION (Q19): What do these results tell us about our three models?
%
% ANSWER (A19): All of our measures indicate that a) the direction
% parameter significantly improves the model and b) model 3 is better than
% both of the other models.

%% Bonus: Kolmogorov-Smirnov plots to evalute models

% Our previous approaches are good for model selection (i.e. adjudicating
% between two models) or for telling us whether a particular covariate
% belongs in the model. However, we would also like to know how well our
% model captures the overall structure of the data. To do this, we will use
% a method related to the Q-Q plot that we explored in an earlier exercise.
% You'll recall that this compares the percentiles of our empirical
% distribution (i.e. the data) to those of some known theoretical
% distribution (e.g. the standard normal). The K-S plot is a similar idea,
% but using the cumulative density functions (CDF) as the basis of
% comparison.

% However, before we do this, we need to transform our spike data so that
% they are identically distributed. That is, we need to somehow take into
% account the fact that lambda is changing over time and as a function of
% our covariates. To do this, we make use of the "time-rescaling theorem,"
% which essentially sums up all of the lambda_t's between each pair of
% spikes to produce rescaled waiting times for each of the spikes. If our
% data are well described by our model, then the re-scaled variables (i.e.
% the Z's in the code below) will be independent, identically distributed
% random variables from the exponential distribution with parameter 1. For
% more details on the math and its application to spiking data see:
%
% Brown EN, Barbieri R, Ventura V, Kass RE, Frank LM. The time-rescaling
% theorem and its application to neural spike train data analysis. Neural
% Comput. 2002 Feb;14(2):325-46. PubMed PMID: 11802915.
% *Available on D2L Course web site

% Count the total # of spikes
nSpikes = length(spikeTimes);

% Evaluate Model 2. This is our model's prediction for the value of lambda
% (our spike rate parameter) in each time bin.
lambda2 = exp(b2(1) + b2(2)*ratPosition + b2(3)*ratPosition.^2);

% Now we rescale waiting times by summing up the lambdas between each pair
% of spikes:
Z = zeros(nSpikes,1);
Z(1) = sum(lambda2(1:spikeIndex(1)));	%1st rescaled waiting time.
for i = 2:nSpikes                       %... and the rest.
  Z(i) = sum(lambda2(spikeIndex(i-1):spikeIndex(i)));
end

% Next we use 'ecdf' to compute the empirical CDF from rescaled waiting times.
[eCDF, zVals] = ecdf(Z);

% Then we compute the theoretical CDF at the same z-values. Recall that the
% theoretical prediction for the rescaled waiting times is an exponential
% distribution with parameter 1. (If you're not sure what this last phrase
% means, go to: https://en.wikipedia.org/wiki/Exponential_distribution)
mCDF = 1-exp(-zVals);                       %Model CDF at z values.

% Now plot the results. Just as for our Q-Q Plots, if the data area well
% described by the theoretical distribution, our plot should fall along the
% y = x line.
figure
subplot(2,1,1);
plot(mCDF,eCDF,'LineWidth',2)                             %Create K-S plot.
hold on

% Use the normal approximation to calculate CIs (good for n > 20)
plot([0 1], [0 1]+1.96/sqrt(nSpikes),'k--')	%Upper confidence bound.
plot([0 1], [0 1]-1.96/sqrt(nSpikes),'k--')	%Lower confidence bound.
xlabel('Model CDF')                         %Label the axes.
ylabel('Empirical CDF')
axis([0,1,0,1]);
title('KS plot of rescaled data for model #2');

% TODO: Evaluate model 3.
lambda3 = exp(b3(1) + b3(2)*ratPosition + ...
    b3(3)*ratPosition.^2 + b3(4)*ratDirection);

% Re-scale the waiting times
Z = zeros(nSpikes,1);
Z(1) = sum(lambda3(1:spikeIndex(1)));	%1st rescaled waiting time.
for i=2:nSpikes							%... and the rest.
  Z(i)=sum(lambda3(spikeIndex(i-1):spikeIndex(i)));
end

[eCDF, zVals] = ecdf(Z);                    %Define empirical CDF,
mCDF = 1-exp(-zVals);                       %...and model CDF

subplot(2,1,2);
h1=plot(mCDF,eCDF,'LineWidth',2);               %...to create KS-plot.
hold on
h2=plot([0 1], [0 1]+1.96/sqrt(nSpikes),'k--'); %Upper 95% confidence bound.
plot([0 1], [0 1]-1.96/sqrt(nSpikes),'k--')     %Lower 95% confidence bound.
xlabel('Model CDF'); ylabel('Empirical CDF');
axis([0,1,0,1]);
title('KS plot of rescaled data for model #3');
legend([h1,h2],{'KS Plot','95% CI'},'Location','Northwest');

% QUESTION (Q20): What do you conclude from a comparison of the two
% Kolmogorov-Smirnov plots?
%
% ANSWER (A20): Model #3 does a reasonably good job of capturing the
% structure in our spike data, whereas model #2 does not. Note that this is
% a slightly more powerful statement than what we could do with, e.g. AIC
% or CI's for our beta parameters. These later metrics are helpful for
% deciding whether one model is better than another in a relative sense,
% but they don't tell us how well our model is doing at capturing the
% overall statisitcal structure of our data in an absolute sense. The KS
% plots do exactly this, so they are useful for reassuring us that our
% chosen framework--i.e the Poisson point process model and a maximum
% likelihood approach--is a good one for the type of data we are trying to
% explain.

%% Bonus: Visualize the model's predicted rate on the raw data:

figure('position',[50 50 1200 600]);
yyaxis left
plot(expTime,ratPosition);
hold on
hp=plot(expTime(spikeIndex),ratPosition(spikeIndex),'k.');
set(hp,'MarkerSize',12);
xlabel('Time (seconds)');
ylabel('Position (cm)')
title('Rat position and predicted spike rate vs. time');
yyaxis right
plot(expTime,lambda3.*1000);
ylabel('Predicted rate (spikes/sec)')
legend(hp,'Spikes','Location','Northwest');

%% Bonus: Form the fit using logistic regression

[b4,dev4,stats4] = glmfit([ratPosition,ratPosition.^2,ratDirection],...
    spikeTrain,'binomial','link','logit');

[pSpike, pCIlow, pCIhi] = glmval(b4,[ratPosition,ratPosition.^2,ratDirection],...
    'logit',stats4);

figure('position',[50 50 1200 600]);
yyaxis left
plot(expTime,ratPosition);
hold on
hp=plot(expTime(spikeIndex),ratPosition(spikeIndex),'k.');
set(hp,'MarkerSize',12);
xlabel('Time (seconds)');
ylabel('Position (cm)')
title('Rat position and predicted spike probability vs. time');
yyaxis right
%plot(expTime,[pSpike,pSpike-pCIlow,pSpike+pCIhi]);
plot(expTime,pSpike);
ylabel('Spiking probability')
legend(hp,'Spikes','Location','Northwest');