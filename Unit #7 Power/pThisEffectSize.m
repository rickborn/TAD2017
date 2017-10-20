function [pVal] = pThisEffectSize(n,dPrime,nSims)

% pThisEffectSize.m: How likely are we to get a given effect size under H0?
%
% Inspired by a failure-to-listen of Slater Sharp

% Assume H0 is true. Draw 'n' samples from same distribution and measure
% dPrime:

allDprime = zeros(nSims,1);
for k = 1:nSims
    allDprime(k) = abs(mean(randn(n,1)) - mean(randn(n,1)));
end

% faster way?
% allSamples = randn(n,2,nSims);
% allDprime = abs(diff(mean(allSamples)));

pVal = sum(allDprime >= dPrime) / nSims;