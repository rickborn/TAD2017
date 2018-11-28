function [D,pHgH,pHgM] = GVTsim(nShots,kStreak,pHit,nSims,pFlag)

% Simulation of bias in counting streaks after Miller & Sanjurjo 2016
% "Surprised by the Gambler's and Hot Hand Fallacies? A Truth in the Law of
% Small Numbers"
%
% D = GVTsim(100,3,0.5,10000,1);
%
% Inputs:
% - nShots: total # of tries in each simulation (default = 100)
% - kStreak: # of hits that constitute a "streak" (default = 3)
% - pHit: probability of a hit on each try (default = 0.5)
% - nSims: # of simulations to run (default = 10,000)
% - pFlag: 1 = plot histogram (default)
%
% Outputs:
% - D: Prob(hit|k hits) - Prob(hit | k misses)
% - pHgH: Prob(hit | k hits)
% - pHgM: Prob(hit | k misses)
%
% Based on counting method and calculations used by:
% Gilovich, T., R. Vallone, and A. Tversky (1985) "The Hot Hand in
% Basketball: On the Misperception of Random Sequences," Cognitive
% Psychology, 17: 295-314.
%
% This is a really nice simulation of a subtle counting bias that went
% undetected for >30 years.
%
% RTB wrote it, 09 November 2018, airplane trip from San Diego to Boston

if nargin < 5, pFlag = 1; end
if nargin < 4, nSims = 10000; end
if nargin < 3, pHit = 0.5; end
if nargin < 2, kStreak = 3; end
if nargin < 1, nShots = 100; end

% init variables to hold results of simulations:
D = zeros(nSims,1);
pHgH = zeros(nSims,1);
pHgH = zeros(nSims,1);
% filter for use with 'conv' to detect streaks of length kStreak
u = ones(1,kStreak);

for k = 1:nSims
    % generate a random sequence of hits/misses:
    x = rand(1,nShots) <= pHit;
    
    % calculate p(hit|k hits)
    % We use the 'conv' trick to find the end of runs of >= kStreak
    w = conv(x,u);
    t = find(w >= kStreak);
    if max(t) == length(x)
        t = t(1:end-1);
    end
    pHitGivenKhits = sum(x(t+1)) / length(t);
    pHgH(k) = pHitGivenKhits;
    
    % calculate p(hit|k misses)
    % same but first convert misses to hits to find miss streaks
    w = conv(~x,u);
    t = find(w >= kStreak);
    if max(t) == length(x)
        t = t(1:end-1);
    end
    pHitGivenKmisses = sum(x(t+1)) / length(t);
    pHgM(k) = pHitGivenKmisses;
    
    % Gilovich et al. calculated the difference, which should be 0 under
    % the null hypothesis:
    D(k) = pHitGivenKhits - pHitGivenKmisses;
end

if pFlag
    histogram(D);
    ylabel('#');
    xlabel(['P(Hit|',num2str(kStreak),' hits) - P(Hit|',num2str(kStreak),' misses)']);
    set(gca,'Fontsize',14);
    % draw vertical lines for the mean and the median:
    Dmean = mean(D,'omitnan');
    Dmedian = median(D,'omitnan');
    ax = axis;
    h1 = line([Dmean,Dmean],[ax(3),ax(4)],'Color',[0.7,0.2,0]);
    h2 = line([Dmedian,Dmedian],[ax(3),ax(4)],'Color',[0,0.7,0.2]);
    legend([h1,h2],{'mean','median'});
    
end
    