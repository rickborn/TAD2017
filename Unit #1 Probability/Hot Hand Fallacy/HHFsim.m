% HHFsim.m: simulation of hot hand fallacy
%
% RTB wrote it, 28 May 2020, North Hero, VT (Lake Champlain)

% Inspired by this Data Colada post by Uri Simonsohn:
% http://datacolada.org/88
%
% My favorite fact for this post: calibration of massiveness of hot hand fallacy
% It is worth noting that expecting r=.40 for consecutive shots
% (see chart above) implies implausibly strong hand-hotness. As a
% calibration, I simulated a player that when cold (previous shot was
% missed) makes as few shots as the NBA player with the lowest shooting
% percentage in the 2019-2020 season (Taurean Prince: 37.6%), but then gets
% dramatically better with each basket that he makes until, after three
% consecutive hits, he acquires the best shooting percentage in the league,
% 74.2%, that of Mitchell Robinson. In a series of 100 shots, this
% impossibly hot-handed player would on average have a correlation of r=.17
% for two consecutive shots

%% First let's calculate the simulated expectation under the null

% We are looking at the degree of correlation for consecutive shots. That
% is, for a sequence of 100 shots, the overall correlation between the two
% vectors of 0s and 1s corresponding to shots 1:99 and 2:100, or, in MATLAB
% speak, if our vector is X: corr([X(1:end-1),X(2:end)])

nSims = 10000;
nShots = 100;
pHit = 0.5;

allR = zeros(nSims,1);
for k = 1:nSims
    x = rand(nShots,1) <= pHit;
    r = corr([x(1:end-1,:),x(2:end,:)]);
    allR(k,1) = r(2,1);
end

figure
histogram(allR);
xlabel('Correlation between consecutive shots');
ylabel('# of simulations');
title('Null Distribution for Hot Hand Fallacy');
nullMean = mean(allR);          % should be around 0

%% Now generate data according to Simonsohn's rule above

nSims = 100000;
nShots = 100;
% pHit = linspace(0.376,0.742,4);
pHit = [0.376, 0.450, 0.508, 0.742];    % vals used by Simonsohn
currStreak = 0;     % keep track of current 'hit' streak

allR = zeros(nSims,1);
for k = 1:nSims
    
    x = zeros(nShots,1);
    for m = 1:nShots
        thisShot = rand < pHit(currStreak+1);
        if thisShot
            x(m,1) = 1;
            currStreak = currStreak + 1;
            if currStreak >= 3
                currStreak = 3; % shooting percentage max's out at 3
            end
        else
            x(m,1) = 0;
            currStreak = 0;
        end
        
    end
    
    r = corr([x(1:end-1,:),x(2:end,:)]);
    allR(k,1) = r(2,1);
end

HHmean = mean(allR);

figure
histogram(allR);
xlabel('Correlation between consecutive shots');
ylabel('# of simulations');
title('Distribution for ''Extreme Hot Handedness''');
hold on

% If I use Simonsohn's values for pHit, I get the same answer of a mean
% correlation of 0.17. If I use a linear spacing from worst to best, I get
% a mean corr of about 0.22.

%% Calculate the 95% CI for our measure

myAlpha = 0.05;

sortedR = sort(allR);
[idxLo,idxHi] = getCIidx(myAlpha,nSims);
ciLo = sortedR(idxLo);
ciHi = sortedR(idxHi);

ax = axis;
xT = ceil(ax(1) * 10) / 10;
yT = ax(4) - (0.1 * ax(4));
text(xT,yT,['Mean corr = ', num2str(HHmean,3)]);

h0 = line([HHmean, HHmean],[ax(3), ax(4)],'Color','y');
h1 = line([ciLo, ciLo],[ax(3), ax(4)],'Color','r');
h2 = line([ciHi, ciHi],[ax(3), ax(4)],'Color','r');
legend(h1,'95% CI');
% Note that the uppper 95% CI is around 0.41, which is the overestimated
% average in the GVT study #4.
