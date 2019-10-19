function [pwr] = plotmany(nPlots,dPrime,n)

if nargin < 1, nPlots = 4; end
if nargin < 2, dPrime = 0.6; end
if nargin < 3, n = 22; end

figure
for k = 1:nPlots^2
    subplot(nPlots,nPlots,k);
    [~,~,pwr] = citest(n,dPrime,0.05,1);
end