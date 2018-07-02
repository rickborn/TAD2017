% STNdemo.m
% 
% This exercise has been adapted from:
% Chapter 10, of "Case studies in neural data analysis" by Kramer & Eden 2016

% Original source of exercise:
% Mark A. Kramer and Uri T. Eden, "Case Studies in Neural Data Analysis",
% MIT Press, 2016, Chapter 10.
% https://mitpress.mit.edu/books/case-studies-neural-data-analysis
% http://github.com/Mark-Kramer/Case-Studies-Kramer-Eden
%
% Adapted by RTB, February 2018

% % Concepts covered:
% 1. Working with rhythmic spike data
% 2. ISI histograms
% 3. Autocorrelation plots
% 4. Spectral estimators for point processes
% 5. History-dependent GLM
% 6. Basis functions for smooth co-variates
% 7. Model selection using AIC, KS plots, analysis of residuals
%
% Recording of a subthalamic nucleus (STN) neuron from a human performing a
% hand movement task
%
% Data: The data consist of single unit spiking activity from one STN
% neuron recorded over 50 trials of a hand movement task. In each trial, a
% patient held a joystick and watched cues appear on a screen. The first
% cue indicated whether the patient should move to the right (1) or left
% (0). The second cue, the GO cue, indicated the time to start moving the
% joystick. The dataset contains two seconds of activity for each trial,
% one second before the GO cue and one second after. We label the period
% before the GO cue comes on as the "planning period," and the period after
% the GO cue the "movement period."
%
% moveDir: a [50 x 1] logical vector of indicators for the movement direction
%          on each of the 50 trials; 1 = right; 0 = left
% trialTime: time axis in msec for each trial; 0 is time of GO cue
% spikeTrain: presence (1) or absence (0) of spike in each 1-ms time bin

%% Load and plot the data
load('STNdata.mat');
% original file name: 'Ch10-spikes-1.mat'

figure('Name', 'Raw data');
imagesc(trialTime,1:50,spikeTrain)	%Construct a rastergram.
colormap(1-gray);					%...and change the colormap.
xlabel('Time (ms)');
ylabel('Trial #');
title('Raster plot of STN cell');

%% Generate separate raster plots for each direction of movement

figure('Name','Raster plots for different directions of movement');
subplot(2,1,1)
imagesc(trialTime,1:25, spikeTrain(~moveDir,:))		%Image left trials,
colormap(1-gray);
xlabel('Time (ms)');
ylabel('Trial #');
title('Left trials');
% Draw a vertical red line at t=0 for GO cue
ax = axis;
line([0,0],[ax(3),ax(4)],'Color','r');

subplot(2,1,2)
imagesc(trialTime,1:25, spikeTrain(moveDir,:))		%...and right trials.
colormap(1-gray);
xlabel('Time (ms)');
ylabel('Trial #');
title('Right trials');
ax = axis;
line([0,0],[ax(3),ax(4)],'Color','r');

%% Plot a peristimulus time histogram

[spikeTimes,~] = find(spikeTrain');         %Find when spikes occur
% Compute histogram. If we want our bins to represent spikes per second per
% trial, we need to divide by the number of trials (50) and convert ms to
% seconds.

% PSTH with 1-ms bins
PSTH = hist(spikeTimes, 1:2000)/size(spikeTrain,1)/0.001;
figure('Name','Peri-stimulus time histogram');
bar(trialTime,PSTH);
xlabel('Time (ms)');
ylabel('Spike rate (spikes/sec)');
title('PSTH with 1-ms bins');

% PSTH with larger bins
msPerBin = 10;
PSTH = hist(spikeTimes, 1:msPerBin:2000);
bar(trialTime(1:msPerBin:2000),(PSTH/size(spikeTrain,1)/0.001)/msPerBin);
xlabel('Time (ms)');
ylabel('Spike rate (spikes/sec)');
tStr = sprintf('PSTH with %d-ms bins', msPerBin);
title(tStr);

%% Compute average spike rates during planning & movement periods

iPlan = find(trialTime < 0);					%Indices for planning.
iMove = find(trialTime >= 0);					%Indices for movement.
%Compute the average spike rate,
% planRate = mean(mean(spikeTrain(:,iPlan)))/0.001;
% moveRate = mean(mean(spikeTrain(:,iMove)))/0.001;

%% Compute average firing rates for left vs. right movements

LRate=mean(mean(spikeTrain(~moveDir,:)))/1e-3;	%... and compute rates.
RRate=mean(mean(spikeTrain(moveDir,:)))/1e-3;

PlanL = sum(spikeTrain(~moveDir,iPlan),2);	%Firing rate L, planning.
PlanR = sum(spikeTrain(moveDir,iPlan),2);	%Firing rate R, planning.
MoveL = sum(spikeTrain(~moveDir,iMove),2);	%Firing rate L, movement.
MoveR = sum(spikeTrain(moveDir,iMove),2);	%Firing rate R, movement.
figure('Name','Mean firing rates for different periods and trial types');
boxplot([PlanL PlanR MoveL MoveR], ...	%Display results.
  'labels',{'Plan Left','Plan Right','Move Left','Move Right'});
ylabel('Firing rate (spikes/sec)');

%% ISI histogram for all trials
ISIs=diff(spikeTimes);	%Determine ISIs for all time & trials.
ISIs=ISIs(ISIs>0);%Remove spurious values between trials,
figure('Name', 'Interspike interval histogram for all trials');
hist(ISIs,0:250);		%...and display results.
ax = axis;
axis([0,250,ax(3),ax(4)]);
xlabel('Interspike interval (ms)');
ylabel('Count');

%% ISI histograms separately for the planning and movement periods

%In the planning period,				%...find spikes,
[spikeTimesPlan, ~]=find(spikeTrain(:,iPlan)');
planISIs = diff(spikeTimesPlan);		%...compute ISIs,
planISIs = planISIs(planISIs>0);	%...drop spurious ones,
figure('Name','ISI histograms');
subplot(2,1,1); hist(planISIs,0:250)		%...and plot it,
ax = axis;
axis([0,100,ax(3),ax(4)]);
xlabel('Interspike interval (ms)');
ylabel('Count');
title('Planning period');

%In the movement period,				%...find spikes,
[spikeTimesMove, ~] = find(spikeTrain(:,iMove)');
moveISIs = diff(spikeTimesMove);		%...compute ISIs,
moveISIs = moveISIs(moveISIs>0);	%...drop spurious ones,
subplot(2,1,2); hist(moveISIs,0:250)		%...and plot it,
ax = axis;
axis([0,100,ax(3),ax(4)]);
xlabel('Interspike interval (ms)');
ylabel('Count');
title('Movement period');

%% Look at the auto-correlation function

nTrials = size(spikeTrain,1);
nLags = size(spikeTrain,2) - 1;
acfPlan = zeros(nTrials,nLags);
acfMove = acfPlan;

for k=1:nTrials								%For each trial,
 planData = spikeTrain(k,iPlan);				%...get planning data,
 moveData = spikeTrain(k,iMove);				%...get movement data,
 acfPlan(k,:)=xcorr(planData-mean(planData),'coeff');%...compute ACFs,
 acfMove(k,:)=xcorr(moveData-mean(moveData),'coeff');
end										%...and plot results,

figure('Name','Autocorrelation of spike increments for lags up to 100 ms');
% By default, 'xcorr' computes the autocorrelation at all lags.
% The lags extend from -1000 ms to +1000 ms. We are only interested in 
% positive lags from 1 ms to 100 ms.
minLag = find(trialTime==1) - 1;
maxLag = find(trialTime == 100) - 1;
subplot(2,1,1); stem(1:100,mean(acfPlan(:,minLag:maxLag)));
ylabel('Autocorrelation')				%... with axes labelled.
title('Planning period');
subplot(2,1,2); stem(1:100,mean(acfMove(:,1001:1100)));
ylabel('Autocorrelation');  xlabel('Lag [ms]')
title('Movement period');

%% Spectral analysis of spiking data

% NOTE: To run the following code, one must download and install the
% Chronux toolbox from chronux.org

%Set the parameters of the MTM.
% The planning and movement periods for each trial are one second in
% duration, so we can achieve a 4 Hz frequency resolution by setting the
% time-bandwidth product to 4.
TW = 4;						%Choose time-bandwidth product of 4.
nTapers=2*TW-1;				%...which sets the # of tapers.
params.Fs=1000;				%Define sampling frequency,
params.tapers=[TW,nTapers];	%...time-band product,# tapers.
params.fpass=[0 500];		%Define frequency range to examine.
params.trialave=1;			%Perform trial averaging.

%Compute the coherence during planning & movement.
[SPlan,~]=mtspectrumpb(spikeTrain(:,iPlan)',params);
[SMove,f]=mtspectrumpb(spikeTrain(:,iMove)',params);

%Plot the spectra ... with axes labelled.
figure('Name','Trial-averaged power spectra');
subplot(2,1,1);
plot(f,SPlan); 
ylabel('Power [Hz]')
title('Planning');
subplot(2,1,2);
plot(f,SMove);
ylabel('Power [Hz]')
xlabel('Frequency [Hz]')
title('Movement');

%% Moving window spectral analysis

% The analysis above assumes that the mean and autocovariance structure do
% not change over time within a period. To see if this is reasonable, we
% will compute a sliding window spectrum.

%Set the parameters of the MTM.
movingwin = [.5 .05];		%Define window duration & step size,
params.fpass = [0 50];		%...frequency range to examine,
params.tapers = [2 3];		%...time-band product, # tapers.
[S,T,F]=mtspecgrampb(spikeTrain',movingwin,params); %Get spectrogram,
T = T + trialTime(1)/1000;	%Set time axis,

figure('Name','Trial-averaged spectogram of spiking data');
imagesc(T,F,S')				%...and display it,
xlabel('Time [s]')			%...with axes labelled.
ylabel('Frequency [Hz]')
c = colorbar;
c.Label.String = 'Power (a.u.)';

%% Poisson point process model #1: movement as the only predictor

% As with many complicated analyses, the heavy lifting is done with a
% single line of MATLAB code (i.e. 'glmfit' below). But sometimes the
% tricky part is formatting the data properly. In this case, we want to
% have one column vector for each co-variate and one
T0 = length(trialTime);                             % # of time points.

% Create an indicator variable for movement. We want zeros for the first
% half of the trial ('planning') and ones for the second half. But then we
% need to repeat this for each of the fifty trials and then concatenate
% them into one long column vector
iMove = ones(nTrials,1)*[zeros(1,T0/2) ones(1,T0/2)];   % [50 x 2000]
% Reshape movement indicator to [100,000 x 1]. The first 2000 points 
% correspond to the first trial. The next 2000 to the 2nd trial, etc. To
% confirm this:
%sum(iMove(1:1000)); % = 0 (first half of first trial; planning period
%sum(iMove(1001:2000)); % = 1000 (2nd half of first trial; movement)
iMove = reshape((iMove)',nTrials*T0,1);

% Also need to do the same thing to the spike data
y = reshape(spikeTrain',nTrials*T0,1);

% Fit model #1
[b1, dev1, stats1] = glmfit(iMove,y,'poisson');

% How do we interpret the parameters in b1? Remember that we're actually
% fitting log(lambda_t), so the first step is to exponentiate. Beta0
% (b1(1)) tells us the predicted firing rate during the planning phase:
planRate = exp(b1(1)) * 1000;   % 38.96 spikes/sec;

% The second parameter tells us how much this baseline rate is modulated by
% the movement
moveMod = exp(b1(2));           % 1.4107, so rate increases by 41%
moveRate = planRate * moveMod;  % 54.96 spikes per second
% Both of these seem reasonable given our previous analyses

% Is the modulation by movement statistically significant? In this case,
% the null value for our movement modulation is 1, so we can ask whether
% the 95% CI includes 1. It does not:
%Compute CI for firing rate modulation during movement period, 
CI_lower = exp(b1(2) - 1.96*stats1.se(2));	%...lower CI,
CI_upper = exp(b1(2) + 1.96*stats1.se(2));	%...upper CI.

% Or we can assess it directly using the Wald test, which we get for free
% from the 'glmfit':
pVal = stats1.p(2);		%p-value from Wald test.

%% KS plot to test how well our model accounts for the data

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

% Evaluate Model 1. This is our model's prediction for the value of lambda
% (our spike rate parameter) in each time bin.
lambda1=exp(b1(1)+b1(2)*iMove);		%Firing rate of Model 1.

% Re-scale the intervals. 
spikeindex=find(y);					%Find spikes,
N=length(spikeindex);				%...and total # spikes.
Z(1)=sum(lambda1(1:spikeindex(1)));	%1st rescaled waiting time,
for i=2:N,							%...and the rest.
     Z(i)=sum(lambda1(spikeindex(i-1)+1:spikeindex(i)));
end

[eCDF, zvals] = ecdf(Z);			%Empirical CDF at z values.
mCDF = 1-exp(-zvals);				%Model CDF at z values.

figure('Name','KS plot for model 1 based on planning vs. movement periods');
plot(mCDF,eCDF,'k-','LineWidth',2);
hold on
plot([0 1], [0 1]+1.36/sqrt(N),'r--');  %Upper confidence bound:
plot([0 1], [0 1]-1.36/sqrt(N),'r--');	%Lower confidence bound.
xlabel('Model CDF');
ylabel('Empirical CDF');
axis([0,1,0,1]);

%% Refining the model: movement direction.

% To get some insight into what we need to add to the model, we can plot
% the residuals as a function of whether the movement was to the right or
% to the left:

%Create an indicator variable for trial movement direction.
iDir = logical(reshape((moveDir*ones(1,T0))',nTrials*T0,1));
i0 = find(iDir==0);		%Left movement trial.
i1 = find(iDir==1);		%Right movement trial.
R = cumsum(stats1.resid);%Cumulative sum of Model 1 residuals.

figure('Name','Point process residuals for model 1');
plot(i0./T0,R(i0),'b.')      %Plot residuals for L trials,
hold on
plot(i1./T0,R(i1),'g.')		%...and plot residuals for R trials,
hL = legend('Left','Right','Location','southeast');
xlabel('Trial #')			%...label axes.
ylabel('Integrate Point Process Residual')

% If our model is good, then the cumulative residuals should bounce around
% randomly but close to zero. As it is, we see that the residuals
% consistently increase for left movement trials, indicating that our model
% is systematically underestimates the rate for these trials, and vice
% versa for right momvement trials. This makes sense, because we saw in our
% earlier plot (fig. 4) that the cell fires more for left than for right
% trials. So adding an indicator variable for direction should improve our
% fit. Let's see.

%% Fit Model 2, and return estimates and useful statistics.

[b2,dev2,stats2] = glmfit([iMove iDir],y,'poisson');

% Is the iDir parameter statistically significant?
pVal = stats2.p(3);

% Repeat KS plot for model #2:
lambda2=exp(b2(1)+b2(2)*iMove+b2(3)*iDir);	%Evaluate Model 2.
Z(1)=sum(lambda2(1:spikeindex(1)));	%1st rescaled waiting time,
for i=2:N							%... and the rest.
  Z(i)=sum(lambda2(spikeindex(i-1)+1:spikeindex(i)));
end
[eCDF, zvals]=ecdf(Z);				%Empirical CDF.
mCDF = 1-exp(-zvals);				%Model CDF.

figure('Name','KS plot for model 2 including trial period and direction');
plot(mCDF,eCDF,'k-','LineWidth',2);
hold on
plot([0 1], [0 1]+1.36/sqrt(N),'r--')		%Upper confidence bound.
plot([0 1], [0 1]-1.36/sqrt(N),'r--')		%Lower confidence bound.
xlabel('Model CDF')					%Label the axes.
ylabel('Empirical CDF')
axis([0,1,0,1]);

%% Refining the model: spiking history dependence.

% Now we will add predictors to our model related to past spiking history.
% But how far back in time should we look? In our spectral anaysis, we
% found a peak during the planning period at 15-20 Hz. A history
% dependent model that extends to lags of 70 ms should be able to capture
% rhythms above 15 Hz. (i.e. 1000/70 = 14.3);
%
% We are thus saying that in order to predict spiking at the current time,
% we need to observe the spiking history 70 ms in the past:

modelOrd = 70;					%Set the model order.
T1 = T0 - modelOrd;				%Update # time points.

%Redefine observable & predictors to support history dependence.
y=reshape(spikeTrain(:,modelOrd+1:end)', nTrials*T1,1);			%Data.
iDir=reshape((moveDir*ones(1,T1))',nTrials*T1,1);             %Direction.
iMove=ones(nTrials,1)*[zeros(1,T0/2-modelOrd) ones(1,T0/2)];	%Period,
iMove=reshape((iMove)',nTrials*T1,1);						%...reshaped.

% Create the history predictor. For each step back in the past, we grab the
% spikes shifted by the appropriate lag. We're just creating a set of
% shifted, truncated copies of our spike train. So our y variable starts at
% time point 71 and goes to the end; our first history colum starts at time
% point 70 and goes to end-1; our 2nd column starts at time point 69 and
% goes to end-2, etc. Note that we truncate each trial before concatenating
% everything into one long column vector.
xHist = zeros(size(y,1),modelOrd);
for i = 1:modelOrd				%...for each step in past.
    %xHist(:,i) = reshape(spikeTrain(:,modelOrd+1-i:end-i)',nTrials*T1,1);
    
    % Doing this in one step makes it harder to parse. This would be
    % clearer:
    thisHistory = spikeTrain(:,modelOrd+1-i:end-i)';
    xHist(:,i) = thisHistory(:);
end

%Fit Model 3, with history dependence.
[b3,dev3,stats3] = glmfit([iMove iDir xHist],y,'poisson');

% However, model #3 assumes that the history effect is independent of the
% movement period. Yet, in our spectral analysis, we observed the strong
% rhythmic activity mainly in the planning phase. Therefore we can improve
% our model by allowing for different effects of history during the
% planning and the movement periods:

%Fit Model 4, with history dependence in each period.
[b4,dev4,stats4] = glmfit([iMove iDir ...       %Period & direction.
  (~iMove*ones(1,modelOrd)).*xHist ...          %History in planning.
  (iMove*ones(1,modelOrd)).*xHist], y,'poisson');%History in movement.

%% KS Plot for model #4

%Evaluate Model 4:
lambda4 = glmval(b4,[iMove iDir ...	
    ((1-iMove)*ones(1,modelOrd)).*xHist ...
    (iMove*ones(1,modelOrd)).*xHist],'log');

spikeindex=find(y);					%Find spikes,
N=length(spikeindex);				%...and total # spikes.
Z(1)=sum(lambda4(1:spikeindex(1)));	%1st rescaled waiting time,
for i=2:N							%... and the rest.
  Z(i)=sum(lambda4(spikeindex(i-1)+1:spikeindex(i)));
end
[eCDF, zvals]=ecdf(Z);				%Empirical CDF.
mCDF = 1-exp(-zvals);				%Model CDF.

figure('Name','KS Plot for model #4');
plot(mCDF,eCDF,'k-','LineWidth',2);	%Create KS-plot.
hold on								%Freeze graphics.
plot([0 1], [0 1]+1.36/sqrt(N),'r--')		%Upper confidence bound.
plot([0 1], [0 1]-1.36/sqrt(N),'r--')		%Lower confidence bound.
xlabel('Model CDF')					%Label the axes.
ylabel('Empirical CDF')
axis([0,1,0,1]);

%% Interpreting our fit parameters

% Model 4 has 143 parameters: one intercept parameter, one parameter for
% the period (planning vs. movement), one for the direction of movement
% (right vs. left), 70 parameters for the history dependence during the
% planning period, and 70 more for the history dependence during the
% momvement period. Whenever we construct a point process model that
% includes spiking history, we call it a "conditional intensity model"
% rather than a rate model.

% We can examine the first three exponentiated parameters of Model 4.
% exp(b4(1:3)) % = [0.0480  1.3819  0.6063]

% Since this model has 140 additional parameters related to spiking history
% dependence, inspecting the values directly would be overwhelming. So
% let's plot their exponentiated values as a function of lag for each
% period:
figure('Name', 'Exponentiated GLM parameters related to history dependence');
subplot(2,1,1); plot(1:modelOrd,exp(b4(4:modelOrd+3)),'LineWidth',2);
hold on
ax = axis;
line([ax(1),ax(2)],[1,1],'Color','r','LineStyle','--');
ylabel('Modulation')							%..axes labelled.
title('Planning');

subplot(2,1,2); plot(1:modelOrd,exp(b4(modelOrd+4:end)),'LineWidth',2);
hold on
ax = axis;
line([ax(1),ax(2)],[1,1],'Color','r','LineStyle','--');
ylabel('Modulation'); xlabel('Lag [ms]')		%..axes labelled.
title('Movement');

% We interpret each exponentiated parameter as a modulation of the spike
% rate when a previous spike has occurred at a given lag. So values greater
% than one indicate an increased probability of spiking, and values less
% than one a decreased probability. For the movement period, we see the
% neuron's refractory period (< 1 for short lags), a peak at 6 ms when
% there is an increased probability of firing, but after 10 ms it hovers
% around one (no modulation). For the planning period, however, we see the
% same early stuff, but then we see the later downward rate modulation
% (15-35 ms after a spike) followed by an upward modulation (45-65 ms)
% corresponding to the band of power between 15-20 Hz in the spectogram.

%% Are these history effects statistically significant?

% To better visualize the p-values associated with each parameter, we plot
% the negative of the natural logarithm of the p-value. Any value above 3
% indicates a p-value < 0.05.

figure('Name','Negative log p-values for model 4 parameters');
subplot(2,1,1)
hP = plot(1:modelOrd,-log(stats4.p(4:3+modelOrd)),'.');	%p-values for Planning,
set(hP,'MarkerSize',15);
hold on
plot([1 modelOrd],[-log(.05) -log(.05)],'r--');     %...draw threshold at p=0.05
ylabel('-Log p-value');
title('Planning');

subplot(2,1,2)					
hP = plot(1:modelOrd,-log(stats4.p(4+modelOrd:end)),'.'); %p-values for Movement
set(hP,'MarkerSize',15);
hold on									%...freeze graphics,
plot([1 modelOrd],[-log(.05) -log(.05)],'r--');	%...draw threshold at p=0.05
ylabel('-Log p-value');
xlabel('Lag (ms)');
title('Movement');

%Compare the two nested GLMs:
pVal = 1-chi2cdf(dev3-dev4,modelOrd);

%% Chapter 10, Choice of model order.

% For the data in the planning period alone, we will find the model order
% that minimizes the Akaike Information Criterion (AIC)
%
% NOTE: This block of code takes several minutes to run!
skipFlag = 1;

if ~skipFlag
    maxOrd = 100;       	%Maximum model order
    
    %Redefine observable & predictors to support history dependence.
    yPlan = reshape(spikeTrain(:,maxOrd+1:T0/2)',nTrials*(T0/2-maxOrd),1);
    planDir = reshape((moveDir*ones(1,T0/2-maxOrd))',nTrials*(T0/2-maxOrd),1);
    planHist=[];			%Create the history predictor,
    aic = zeros(maxOrd,1);
    
    for i=1:maxOrd			%...for each step in past,
        planHist = [planHist ...	%...define history, concatenate for each higher order
            reshape(spikeTrain(:,maxOrd+1-i:T0/2-i)',nTrials*(T0/2-maxOrd),1)];
        %..fit the model,
        [b0,dev0,stats0]=glmfit([planDir planHist],yPlan,'poisson');
        %...and compute the AIC.
        aic(i) = dev0 + (2*length(b0));
    end
    
    figure('Name',...
        'AIC plot of models with increasing number of history-dependent components');
    plot(1:maxOrd,aic,'b-','LineWidth',2);     %Plot the AIC,
    xlabel('Model Order')   %...with axes labelled,
    ylabel('AIC')
    
    allOrd = 1:maxOrd;
    bestOrder = allOrd(aic == min(aic));    % 'best' = model with 62 ms history
    
end

%% Refining the model: smooth history dependence.

% Do we really believe that the influence of a spike 25 ms in the past can
% can be very different from one 26 ms in the past? If not, it might be
% better to define a model with fewer free parameters that guarantees a
% smooth modulation as a function of lag. To do this, we replace the
% separate terms for each single millisecond lag with a small set of basis
% functions on the spiking history.

% Construct Gaussian kernels to use as basis functions:
X = -5:10:modelOrd;
C = zeros(modelOrd,length(X));
for i=1:modelOrd		
    C(i,:) = normpdf(X,i,5);
end
% To visualize the kernels: plot(C)

% The result in C is a [70 x 8] matrix that converts the previous 70-column
% xHist matrix into a new seet of 8 columns to include in the design
% matrix: xHist*C. Cute! So now, instead of 70 free parameters to capture
% spike history, we reduce it to only 8.

%Fit Model 5, with Gaussian kernel basis.
nParams = size(C,2);

[b5,dev5,stats5] = glmfit([iMove iDir ...	%Period & direction.
  (~iMove*ones(1,nParams)).*(xHist*C) ...   %History in planning.
  (iMove*ones(1,nParams)).*(xHist*C)], ...	%History in movement.
  y,'poisson');

% Now let's replot the fit parameters that correspond to our history terms:
figure('Name','History fit parameters using Gaussian kernel method');
subplot(2,1,1);
plot(1:modelOrd,exp(C*b5(4:nParams+3)),'b-','LineWidth',2);    %Planning
hold on
ax = axis;
line([ax(1),ax(2)],[1,1],'Color','r','LineStyle','--');
ylabel('Modulation')
title('Planning');

subplot(2,1,2);
plot(1:modelOrd,exp(C*b5(nParams+4:end)),'b-','LineWidth',2);    %Movement
hold on
ax = axis;
line([ax(1),ax(2)],[1,1],'Color','r','LineStyle','--');
ylabel('Modulation'); xlabel('Lag (ms)')
title('Movement');

%% KS Plot for model 5

%Evaluate Model 5.
lambda5 = glmval(b5,[iMove iDir ...	
  (~iMove*ones(1,nParams)).*(xHist*C) ...       % planning phase
  (iMove*ones(1,nParams)).*(xHist*C)],'log');   % movement phase

% Re-scale intervals
spikeindex=find(y);					%Find spikes,
N=length(spikeindex);				%...and total # spikes.
Z(1)=sum(lambda5(1:spikeindex(1)));	%1st rescaled waiting time,
for i=2:N							%... and the rest.
  Z(i)=sum(lambda5(spikeindex(i-1)+1:spikeindex(i)));
end;
[eCDF, zvals]=ecdf(Z);				%Empirical CDF.
mCDF = 1-exp(-zvals);				%Model CDF.

figure('Name','KS Plot for model 5: history modulation with Gaussian kernels');
plot(mCDF,eCDF,'k-','LineWidth',2);	%Create KS-plot.
hold on
plot([0 1], [0 1]+1.36/sqrt(N),'r--')		%Upper confidence bound.
plot([0 1], [0 1]-1.36/sqrt(N),'r--')		%Lower confidence bound.
xlabel('Model CDF')					%Label the axes.
ylabel('Empirical CDF')
axis([0,1,0,1]);

%% Evaluating model #5 with a plot of p-values

figure('Name','Negative log p-values for model 5 history parameters');

% Correct our p-value threshold for multiple comparisons:
myAlpha = 0.05;
pThresh = myAlpha / nParams;

subplot(2,1,1)							%p-values for Planning,
hP = plot(1:nParams,-log(stats5.p(4:3+nParams)),'.');
set(hP,'MarkerSize',15);
hold on									%...freeze graphics,
plot([1 nParams],[-log(pThresh) -log(pThresh)],'r--')	%...draw threshold,
ylabel('-Log p-value');
title('Planning');

subplot(2,1,2)							%p-values for Movement,
hP = plot(1:nParams,-log(stats5.p(4+nParams:end)),'.');
set(hP,'MarkerSize',15);
hold on									%...freeze graphics,
plot([1 nParams],[-log(pThresh) -log(pThresh)],'r--')	%...draw threshold,
ylabel('-Log p-value');
xlabel('Lag parameter');
title('Movement');

%% Compare different models for planning vs movement phase

% Fit a reduced version of Model 5, where spike history and period are
% independent. That is, one where we assume that the effect of spiking
% history is the same in both phases.
[b6,dev6,stats6] = glmfit([iMove iDir xHist*C],y,'poisson');

pVal = 1-chi2cdf(dev6-dev5,nParams);	%Compare Model 5 and Model 6.

% This tells us that model 5 really is much better, and we can be confident
% that the cell exhibits different spiking dynamics during the two phases
% of the task (planning vs. movement). The rhythmic structure that we see
% during planning is attenuated during the movement period.