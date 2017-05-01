function [p] = ttest_perm(A,B,nSims)

% Exercise: Place comments at the '%' in the following lines of code.
% The general goal of the code is to perform a permutation version of a
% 2-sample t-test. The two samples are contained in 'A' and 'B'. 'nSims' is
% the number of simulations to run in order to get a p-value.

% RTB wrote it, 19 Aug 2012

% setting a default value for # of simulations
if nargin < 3, nSims = 10000; end

% 
C = [A(:); B(:)];

% 
realDiff = mean(B) - mean(A);

% 
bsDiffs = zeros(1,nSims);

% 
for iSim = 1:nSims
    bsC = C(randperm(numel(C)));    % 
    bsA = bsC(1:numel(A));          % 
    bsB = bsC(numel(A)+1:end);      % 
    bsDiffs(iSim) = mean(bsB) - mean(bsA);  % 
end

% 
if realDiff <= mean(bsDiffs)
    t = find(bsDiffs <= realDiff);
else
    t = find(bsDiffs >= realDiff);
end

% 
p = length(t) / nSims;