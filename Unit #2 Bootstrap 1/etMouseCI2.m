% etMouseCI2.m
%
% Bootstrap confidence intervals for one-sample
% Ch. 13 of E & T
%
% RTB wrote it, 04 Jan. 2017

% Concepts covered:
% 1. Confidence intervals from standard normal distribution
% 2. CIs by bootstrap: percentile method
% 3. Legends using plot handles

% from table 2.1
rxGrp = [94,197,16,38,99,141,23];   % mean of 86.857, sigma 66.767
nRx = length(rxGrp);

myAlpha = 0.10;     % E&T use 90% CIs in their examples

muRx = mean(rxGrp);
seRx = std(rxGrp) / sqrt(nRx);

% to calculate the standard confidence interval, we just use the normal
% distribution
CI90lo = muRx + (seRx*(norminv(myAlpha/2)));
CI90hi = muRx + (seRx*(norminv(1 - myAlpha/2)));

% Bootstrap the means
nBoot = 100000;
muBS = zeros(nBoot,1);
% NOTE that we are not just looking at the distribution of our means, but
% rather caculating a bootstrap z-statistic for each bootstrap sample.
% See p. 160 of E & T
for k = 1:nBoot
    xStar = rxGrp(unidrnd(nRx,nRx,1));
    muBS(k) = mean(xStar);    
end
seRxBS = std(muBS) * sqrt(nRx/(nRx-1)); % see p. 14 of E&T for 'correction factor'

% Teachable moment: the BS derivation of std. error doesn't take account of
% the unbiased estimate of variance, which divides the summed squared
% errors by (n-1) instead of n. You can see this by computing the biased
% est. of var using 'var(x,1)' (forces normalization by n instead of n-1)
% sqrt(var(rxGrp,1)/nRx); ans = 23.3635
% std(muBS); ans = 23.3927

% Use built-in MATLAB function 'bootci' to do the same thing
[CI90] = bootci(nBoot,{@mean,rxGrp},'alpha',myAlpha,'type','percentile');

% Fig. 13.1 distribution of our bootstrapped values:
figure, hist(muBS,100);
xlabel('Bootstrap estimate of mean'); ylabel('#');
title('Figure 13.1 of E&T, p. 169');

sortedMUbs = sort(muBS);
idxLo = nBoot * myAlpha/2;
idxHi = nBoot * (1 - myAlpha/2);
CI90BSlo = sortedMUbs(idxLo);
CI90BShi = sortedMUbs(idxHi);

% plot bootstrap CIs on histogram of means
ax = axis;
hL1 = line([CI90BSlo,CI90BSlo],[ax(3),ax(4)],'Color','r');
line([CI90BShi,CI90BShi],[ax(3),ax(4)],'Color','r');

% plot standard normal CIs
hL2 = line([CI90lo,CI90lo],[ax(3),ax(4)],'Color','g');
line([CI90hi,CI90hi],[ax(3),ax(4)],'Color','g');

% plot the mean
hL3 = line([muRx,muRx],[ax(3),ax(4)],'Color','y');

legend([hL1,hL2,hL3],'90% CI by Bootstrap','90% CI by Standard Normal',...
    'Mean','Location','northeast');