function [y,y1D] = JSestimator(x,N)

% JSestimator.m: The James-Stein empirical Bayes estimator, from Efron
%
% y = muJS(x,N,pFlag)
%
% Inputs:
% - x, a vector of the individual means for a number of samples
% - N, for binomial data, the number of observations that went into each mean
%
% Outputs:
% - y, the James-Stein empirical Bayes estimator
%
% see pp. 6-7 of Efron 2010
%
% RTB wrote it, 10 Feb. 2014 (plane to SPC Study Section in San Francisco)

if nargin < 2
    varX = var(x);
else
    varX = mean(x).*(1-mean(x)) ./ N;   % uses binomial estimate of variance
end

% The wgtFactor determines how much we pull our estimate towards the Grand
% Mean (wgtFactor = 0) vs. towards the MLE (individual means; WF=1). It is
% telling us how much weight we should give to the residuals.
wgtFactor = 1 - (((length(x)-3) .* varX) / sum((x - mean(x)).^2));
y = mean(x) + (wgtFactor .* (x - mean(x)));

