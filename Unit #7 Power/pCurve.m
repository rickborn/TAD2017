function [statPower,perCentPerBin,xVals] = pCurve(dPrime,n,nSim,colStr)

% pCurve.m: generates p-curves as proposed by Simonsohn et al. 2014
%
% e.g. [statPower,perCentPerBin,xVals] = pCurve(0.57,50,10000,'b');
%
% Inputs:
% - dPrime, a measure of effect size: d' = (M2 - M1) / SD (default = 0.64)
% - n, sample size (default = 20)
% - nSim, number of simulations to run (default = 10,000)
% - colStr, color to use for plot (default = red)
%
% Outputs:
% - statPower, statistical power (1 - beta)
% - perCentPerBin, y-axis of graph: %age of p-values at a given criterion
% - xVals, centers of bins for the plot
%
% Simonsohn U, Nelson LD, Simmons JP. p-Curve and Effect Size: Correcting for
% Publication Bias Using Only Significant Results. Perspect Psychol Sci. 2014
% see fig. 1, p. 668
%
% see also: pCurveFamily.m
%
% RTB wrote it, May 2015

if nargin < 4, colStr = 'r'; end
if nargin < 3, nSim = 10000; end
if nargin < 2, n = 20; end
if nargin < 1, dPrime = 0.64; end

rng('shuffle');

% d-prime is measure of effect size: d' = (M2 - M1) / SD
%
% for 2-sample t-test: t = (M2-M1) / (SD * sqrt(2/n))
% d' = t * sqrt(2/n)

% Generate samples under HA for the specificed value of d-prime
S1 = normrnd(1,1,n,nSim);
S2 = normrnd(1+dPrime,1,n,nSim);

% Do all of our t-tests at once:
[~,pVals] = ttest2(S2,S1);

% Idea is to base test on *published* p-values, so we are interested in
% those <= 0.05. We want to know what percentage of p-values fall within
% intervals of:
% 
% [0 to 0.01],[0.01 to 0.02],[0.02 to 0.03], [0.03 to 0.04], [0.04 to 0.05]
%
% To follow the convention of p-Curve, we will plot the percentage of
% values within each interval above the upper edge of that interval:
xVals = [0.01:0.01:0.05];
binCtrs = xVals - 0.005;
nPerBin = hist(pVals(pVals < 0.05),binCtrs);
perCentPerBin = (nPerBin ./ sum(nPerBin)) .* 100;

% Calculate the statistical power of our test.
% Power is 1 - P(type 2 error)
% Since H0 is by definition false, all pVals > 0.05 are type 2 errors.
statPower = 1 - (sum(pVals >= 0.05) / nSim);

% Plot the results
plot(xVals,perCentPerBin,[colStr,'-']); hold on;
h = plot(xVals,perCentPerBin,[colStr,'o']);
set(h,'MarkerFaceColor',colStr);
ax = axis;
axis([0.01 0.05 0 80]);
xlabel('p-value'); ylabel('% of p-values');

% Label each point with its value:
yText = num2str(round(perCentPerBin'));
text(xVals,perCentPerBin+3,yText);

% Print information about this simulation:
gStr = sprintf('n = %d\nd-prime = %0.2f\npower = %d%%\n# sims = %d',...
                n,dPrime,round(statPower*100),nSim);
text(0.025,60,gStr);
