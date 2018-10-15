function [pwr] = plotmany(nPlots,dPrime)

if nargin < 1, nPlots = 3; end
if nargin < 2, dPrime = 0.6; end

figure
for k = 1:nPlots^2
    subplot(nPlots,nPlots,k);
    [~,~,pwr] = citest(22,dPrime,0.05,1);
end