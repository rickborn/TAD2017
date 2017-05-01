function [statPower,perCentPerBin,xVals] = pCurve(dPrime,n,nSim)

% pCurve.m: generates p-curves as proposed by Simonsohn et a. 2014
%
% Inputs:
% - dPrime, a measure of effect size: d' = (M2 - M1) / SD (default = 0.64)
% - n, sample size
% - nSim, number of simulations to run
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

[~,pVals] = ttest2(S2,S1);

% Idea is to base test on *published* p-values, so we are interested in
% those <= 0.05
xVals = [0.00:0.01:0.05];

for iVal = 1:length(xVals)-1
    nPerBin(iVal) = sum(pVals > xVals(iVal) & pVals < xVals(iVal+1));
end
perCentPerBin = (nPerBin ./ sum(nPerBin)) .* 100;
% h = bar(xVals,perCentPerBin);
xVals = xVals(2:end);
plot(xVals,perCentPerBin,'k-'); hold on;
h = plot(xVals,perCentPerBin,'ko');
set(h,'MarkerFaceColor','k');
ax = axis;
axis([0.01 0.05 0 80]);
xlabel('p-value'); ylabel('Percent of p-values');

yText = num2str(round(perCentPerBin'));
text(xVals,perCentPerBin+3,yText);

% Power is 1 - P(type 2 error)
% Since H0 is by definition false, all pVals > 0.05 are type 2 errors
statPower = 1 - (sum(pVals >= 0.05) / nSim);

%tStr = sprintf('n = %d, d-prime = %0.2f, power = %d%%', n,dPrime,round(statPower*100));
%title(tStr);

gStr = sprintf('n = %d\nd-prime = %0.2f\npower = %d%%', n,dPrime,round(statPower*100));
text(0.025,60,gStr);
