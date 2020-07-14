% NDflood.m
%
% Quick probability calculation from Nate Silver's book
%
% RTB wrote it, 14 July 2019, Kinsale, Ireland; end of Cork Distance Week

% Taken from Chapter 6 of Nate Silver's "The Signal and the Noise"
%
% Heading: "The Importance of Communicating Uncertainty"
%
% "In April 1997, the Red River of the North flooded Grand Forks, North
% Dakota, overtopping the town's levees and spilling more than two miles
% into the city. Although there was no loss of life, nearly all of the
% city's 50,000 residents had to be evacuated, cleanup costs ran into the
% billions of dollars, and 75 percent of the city's homes were damaged or
% destroyed.
%
% . . . Residents of Grand Forks had been aware of the flood threat for
% months. Snowfall had been especially heavy in the Great Plains that
% winter, and the National Weather Service, anticipating runoff as the snow
% melted, had predicted the waters of the Red River would crest to
% forty-nine feet, close to the all-time record.
%
% There was just one small problem. The levees in Grand Forks had been
% built to handle a flood of fifty-one feet. Even a small miss in the
% forty-nine-foot prediction could prove catastrophic.
%
% . . . The margin of error on the Weather Service's forecast--based on how
% well their flood forecasts had done in the past--was about plus or minus
% nine feet . . . ."

% QUESTION 1: Assuming what Silver is calling the "margin of error"
% corresponds to a 95% confidence interval, what is the probability that
% the water exceeded the 51-foot height of the levee?

p1 = 1 - normcdf(51,49,9/1.96);
% ANSWER: 0.33
%
% Explanation: In order to use 'normcdf' we need to convert the confidence
% interval to a standard deviation--recall that the standard error of the
% mean is simply the standard deviation of the sampling distribution of the
% mean. You should carry around in your head the so-called "68-95-99.7
% rule," which describes the probabilty mass contained within plus-or-minus
% 1, 2 and 3 standard deviations from the mean, respectively. So a 95%
% confidence interval corresponds to about 2 standard deviates (actually,
% it's 1.96). Also recall that 'normcdf' gives us the area of the normal
% distribution less-than-or-equal-to a given value (in this case, the
% height of the levees, or 51 feet), so to find the probability of the
% water level exceeding 51 feet, we subtract the area from the total area
% of 1. We could also use:
p2 = normcdf(51,49,9/1.96,'upper');

% Suppose you didn't know about 'normcdf' (or, for that matter, 'normpdf',
% but you had extensive experience with 'randn'.
nSim = 100000;
% Normal distribution with a mean of 49 and an sd of 4.6:
mySimHgts = (randn(nSim,1) .* (9/1.96)) + 49;
pSim = sum(mySimHgts > 51) / nSim;

% Follow-up (Silver's book): "In fact, the river crested to fifty-four
% feet."

% QUESTION 2: What is the probability if the "margin of error" corresponds
% to the standard error?

p3 = normcdf(51,49,9,'upper');
% ANSWER: 0.41

% Instead of giving away the error, let's calculate it from some fake
% historical data on the predicted and measured flood heights:
ds = readtable('FloodHgtHistorical.xlsx');

% simulate a bias: NWS adds a 2-foot fudge-factor to their prediction:
%ds.Predicted = ds.Predicted + 2;

% look at the data:
fig1 = figure('Position',[50 10 600 900],'Name','Historical Flood Heights in North Dakota');
subplot(2,1,1);
plot(ds.Year,ds.Measured,'ro');
hold on
xlabel('Year');
ylabel('Flood Hgt (ft)');
plot(ds.Year,ds.Predicted,'b*');
legend('Measured','Predicted');

% look at the distribution of errors:
subplot(2,1,2);
histogram(ds.Predicted - ds.Measured);
xlabel('Prediction Error (ft)');
ylabel('# of Years');

% calculate the measurement error:
mySD = std(ds.Predicted - ds.Measured);
p4 = 1 - normcdf(51,49,mySD);
% ANSWER: 0.336

% BONUS: calculate confidence intervals for the parameters:
myAlpha = 0.05;     % corresponds to 95% CI
[muHat,sigmaHat,muCI,sigmaCI] = normfit(ds.Predicted - ds.Measured);

% QUESTION: Can we reject the null hypothesis that the predictions are not
% biased?

% ANSWER: In this case, a bias would correspond to the mean difference in
% the predicted - measured to be different from 0 (i.e. the predicted tends
% to be significantly greater than or less than the measured). For our
% historical data, we see that the 95% CI for the mean difference includes
% the null value (zero), so we fail to reject H0 and conclude that our data
% do not support the concern that our predictions are systematically
% biased.

% BONUS: 95% CI for the data done two ways:
%
% Method #1: normal assumption (NA):
myZ = norminv(1 - myAlpha/2);
ciNAhi = muHat + (myZ * sigmaHat);
ciNAlo = muHat - (myZ * sigmaHat);

% draw these on our plot
ax = axis;
h1 = line([ciNAlo, ciNAlo],[ax(3),ax(4)],'Color','r');
h2 = line([ciNAhi, ciNAhi],[ax(3),ax(4)],'Color','r');

% Method #2: empirical (EM):
n = length(ds.Predicted);
sortedErr = sort(ds.Predicted - ds.Measured);
idxLo = floor((myAlpha/2) * n);     % index corresponding to lower bound
idxHi = ceil((1-(myAlpha/2)) * n);  % index corresponding to upper bound
ciEMhi = sortedErr(idxHi);
ciEMlo = sortedErr(idxLo);

