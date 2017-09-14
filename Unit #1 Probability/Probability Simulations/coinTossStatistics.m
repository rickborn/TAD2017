function [maxRuns, nTransitions] = coinTossStatistics(nTosses,nSims,pFlag)

% coinTossStatistics: simulates coin tosses and calculates runs and
% switches
%
% [maxRuns, nTransitions] = coinTossStatistics(nTosses,nSims,pFlag)
%
% Inputs:
% - nTosses: # of tosses to simulate in each trial (default = 100)
% - nSims: # of trials to simulate (default = 2000)
% - pFlag: 1 = plot results (default)
%
% Outputs:
% - maxRuns: length of the maximum run in each trial
% - nTransitions: # of switches (H-to-T or vice versa) in each trial
%
% from Gelman & Nolan, Teaching Statistics (A Bag of Tricks), 8.3.2, Real
% vs. fake coin flips. pp. 119-121
%
% Tests students' intuition about random sequences
% Idea is to divide class in half and have one group generate a sequence of
% 0s and 1s by actually tossing a coin 100 times and the other to come up with their
% own random sequence of 0s and 1s without using a coin or a rng.
%
% RTB wrote it, 10 September 2017, about to go searching for Memphre

if nargin < 3, pFlag = 1; end
if nargin < 2, nSims = 2000; end
if nargin < 1, nTosses = 100; end

% random sequences of 0s and 1s
allTosses = round(rand(nTosses,nSims));

allTransitions = abs(diff(allTosses));
nTransitions = sum(allTransitions);

% This works for a vector, but 'find' will not work along colums.
% maxRuns = max(diff(find(allTransitions)));
maxRuns = zeros(1,nSims);
for k = 1:nSims
    maxRuns(k) = max(diff(find(allTransitions(:,k))));
end

if pFlag
    % jitter each variable randomly over [-0.25 to +0.25]
    xJitter = ((rand(size(nTransitions)) - 0.5)/2) + nTransitions;
    yJitter = ((rand(size(maxRuns)) - 0.5)/2) + maxRuns;
    %figure, plot(xJitter,yJitter,'k.');
    figure
    scatterhist(xJitter,yJitter,'Kernel','on','Location','Northwest',...
        'Direction','out','Marker','.');
    hold on
    hp = plot(median(nTransitions),median(maxRuns),'rs','LineWidth',2,'MarkerSize',10);
    xlabel('Number of switches');
    ylabel('Length of longest run');
    tStr = sprintf('%d simulations of %d coin tosses',nSims,nTosses);
    title(tStr);
    legend('Individual sims','Median');
end

end

