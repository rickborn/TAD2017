function [pow] = powSimTtest2(H0,HA,nSims)

% posSimTtest2.m: calculation of power for a 2-sample t-test by simulation
%
% [pow] = powSimTtest2(H0,HA,nPerGroup,nSims)
%
% Inputs:
% - H0: 3-vector containing mean,SD & n for null distribution
% - HA: 2-vector containing mean,SD & n for alternate distribution
% - nSims: number of simulations to run
% 
% RTB wrote it, Biostats cert course, then updated on 24 Sept. 2016
%
% See also: code in 'C:\usr\rick\doc\Committees\PIN\PIN Director\Courses\Stats\Biostatistics 2015-16\MATLAB code\Power and Sample Size'
% See also: sampsizepwr

if nargin < 3, nSims = 100000; end
if nargin < 2, HA = [98,20,64]; end
if nargin < 1, H0 = [85,20,64]; end

allH = zeros(nSims,1);

% For loop method (slow)
% tic; pow = powSimTtest2([85,20,50],[95,20,50],100000); toc
% Elapsed time is 25.151586 seconds.
% for iSim = 1:nSims
%     x = normrnd(H0(1),H0(2),H0(3),1);
%     y = normrnd(HA(1),HA(2),HA(3),1);
%     allH(iSim) = ttest2(x,y);
% end

% vectorized method (fast)
% tic; pow = powSimTtest2([85,20,50],[95,20,50],100000); toc
% Elapsed time is 0.747981 seconds.

allH = ttest2(normrnd(H0(1),H0(2),H0(3),nSims), ...
    normrnd(HA(1),HA(2),HA(3),nSims));

pow = sum(allH) / nSims;

% compare with built-in MATLAB command
% pwr = sampsizepwr('t2',[85 20],95,[],64)
% pwr = 0.8015
%
% or with built-in STATA command:
% power twomeans 85 95, sd1(20) sd2(20) n1(64) alpha(0.05)
% answer: power = 0.8015

% If you want to pick a sample size, you can just plot a curve of the power
% at different values of n.
% allN = [40:5:80];
% allPow = zeros(length(allN),1);
% allN = allN';
% for k = 1:length(allN)
%     allPow(k) = powSimTtest2([85,20,allN(k)],[95,20,allN(k)],100000);
% end
% plot(allN,allPow,'b-');
% xlabel('N'); ylabel('Power');