function [pValCDF,pValPDF] = pNflips2getzHeads(zHeads,N,pHead,nSims,pFlag)

% pNflips2getzHeads.m: probability of taking N or fewer flips to get z heads
%
% What is the probability of taking N or fewer flips to get z heads?
% e.g. [pVal,pValPDF] = pNflips2getzHeads(8,25,0.5,100000,1);
%
% Inputs:
% - zHeads: number of heads to achieve (terminate flipping)
% - N: number of flips we actually needed to get z heads
% - pHead: the probability of a head (default = 0.5)
% - nSims: # of simulations to run (default = 10,000)
% - pFlag: 1 = plot histograms (default = 0)
%
% Outputs:
% - pValCDF: probability of taking N or fewer flips to get z heads
% - pValPDF: probability of taking exactly N flips to get z heads
%
% RTB wrote it, 30 Oct. 2017, after the big wind/rain storm
%
% see also: 'nbincdf' and 'nbinpdf'

% NOTE: This exercise was inspired by Ch. 11 of John Kruschke's "Doing
% Bayesian Data Analysis" (2nd Ed.) in which he points out that how we
% calculate our null model depends on our intent in doing the data
% collection. He gives the example of asking an assistant to flip a coin a
% few times as we watch. The following sequence is obtained:
%
%    TTHHTTHTTTTTTTTTHTTHHTTH
%
% That is, we observe N=24 flips and get z=7 heads. We would like to derive
% the probability of getting a proportion of heads that is 7/24 or smaller
% if the null hypothesis (pHeads=0.5) is true.
%
% Suppose we ask the assistant why he stopped flipping the coin, and he
% says that his lucky number is 24, so he decided to stop when he completed
% 24 flips. In this case, it is easy to calculate a p-value from the
% binomial distribution:
% p = binocdf(7,24,0.5);  % p = 0.032
%
% But suppose instead that the assitant replied that 7 was his lucky
% number, and he stopped when he got 7 heads. This is a very different
% data-generating process and will give us a different p-value. We can
% calculate it by simulation (as below), or we can use the negative
% binomial distribution, y = nbincdf(x,R,p), which gives the probability of
% x-or-fewer failures *before* we get R successes given a single-trial
% probability of success of p. In our case, we want to know the probability
% of 17 or more failures before 7 successes. In MATLAB, this would be:
% p = 1 - nbincdf(16,7,0.5);    % p = 0.0173
%
% or, more accurately (see 'doc nbincdf'):
% p = nbincdf(16,7,0.5,'upper'); % p = 0.0173
%
% y = nbincdf(x,R,p,'upper') returns the complement of the negative
% binomial cdf at each value in x, using an algorithm that more accurately
% computes the extreme upper tail probabilities.

if nargin < 5, pFlag = 0; end
if nargin < 4, nSims = 10000; end
if nargin < 3, pHead = 0.5; end
if nargin < 2, N = 24; end
if nargin < 1, zHeads = 7; end

allProps = zeros(nSims,1);
allN = allProps;
for k = 1:nSims
    nFlips = 0;
    nHeads = 0;
    while nHeads < zHeads
        nFlips = nFlips + 1;
        nHeads = nHeads + (rand < pHead);
    end
    allProps(k) = nHeads / nFlips;
    allN(k) = nFlips;
end
pValCDF = sum(allProps >= (zHeads ./ N)) ./ nSims;
pValPDF = sum(allProps == (zHeads ./ N)) ./ nSims;

if pFlag
    figure, histogram(allProps);
    xlabel('Sample proportion (z/N)');
    ylabel('#');
    tStr = sprintf('Number of flips to get %d heads if p(H) = %.2f',...
        zHeads,pHead);
    title(tStr);
    ax = axis;
    line([zHeads/N,zHeads/N],[ax(3),ax(4)],'Color','r');
    pStr = sprintf('p = %.3f',pValCDF);
    text(ax(2)*0.8,ax(4)*0.9,pStr);
    
    figure, histogram(allN);
    xlabel('# of Flips');
    ylabel('#');
    title(tStr);
    ax = axis;
    line([N,N],[ax(3),ax(4)],'Color','r');
    text(ax(2)*0.8,ax(4)*0.9,pStr);
end