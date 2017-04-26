function [pVal] = FourChoiceProbabilityFunction(nFolks,R,k,nSims,numFlag)

% 4-choice probability: probability of k or more choosing a given #
%
% [pVal] = FourChoiceProbabilityFunction(nFolks,R,k,nSims,numFlag)
%
% For nFolks people "randomly" picking a number from 1 to R (inclusive) 
% what is the probability that k or more folks choose the same number?
%
% Inputs:
% - nFolks: # of people participating in the game (default = 65)
% - R: range of numbers picked from, 1 to R, inclusive (default = 4)
% - k: number of folks picking a given number (default = 28)
% - nSims: number of simulations to run (default = 100,000)
% - numFlag: set for a particular number (default = 0)
%
% Outputs:
% - pVal: probability of k or more folks choosing the same number
%
% Note that this same code can be used to calculate probabilities for the
% birthday problem. E.g. If we have 23 people in a room, what is the
% probability that 2 or more of them share a birthday?
% FourChoiceProbabilityFunction(23,365,2,100000,0)
%
% RTB wrote it

if nargin < 5, numFlag = 0; end
if nargin < 4, nSims = 100000; end
if nargin < 3, k = 28; end
if nargin < 2, R = 4; end
if nargin < 1, nFolks = 65; end
% nFolks = 65; R = 4; nSims = 100000; k = 28;

allSims = unidrnd(R,nFolks,nSims);
allBinned = hist(allSims,[1:R]);

if numFlag
    allSuccesses = allBinned(numFlag,:) >= k;     % only a specific number
else
    allSuccesses = any(allBinned >= k);    % 1, 2, 3 or 4
end

pVal = sum(allSuccesses) / nSims;
