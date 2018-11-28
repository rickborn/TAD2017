function [npAlpha,npBeta] = npSim(crit,dPrime,nSamp,nSims)

% npSim.m: simulates the Neyman-Pearson hypothesis testing algorithm
%
% Inspired by pp.18-20 of Efron's "Computer Age Statistical Inference"
% RTB wrote it, 28 October 2018, in DC for Frank's 60th

% This is currently not giving me the exact same answers as in the Efron
% textbook, but the general shape of the curve is correct.

if nargin < 4, nSims = 100000; end
if nargin < 3, nSamp = 10; end
if nargin < 2, dPrime = 0.5; end

% for tallying false positive and false negatives
fpCount = 0;
fnCount = 0;

for k = 1:nSims
    % draw two samples from two different distributions: one null, one
    % non-null
    S0 = randn(nSamp,1);
    S1 = randn(nSamp,1) + dPrime;
    
    % calculate the likelihood ratio for the null sample
    L0 = sum(log10(normpdf(S0,dPrime,1))) - sum(log10(normpdf(S0,0,1)));
    
    % Type 1 error: choosing f1 when the data came from f0
    if L0 >= crit
        fpCount = fpCount + 1;
    end
    
    % calculate the likelihood ratio for the non-null sample
    L1 = sum(log10(normpdf(S1,dPrime,1))) - sum(log10(normpdf(S1,0,1)));
    
    % Type 2 error: choosing f0 when the data came from f1
    if L1 < crit
        fnCount = fnCount + 1;
    end
end

npAlpha = fpCount / nSims;
npBeta = fnCount / nSims;

end

