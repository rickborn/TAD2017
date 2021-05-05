% IQtest.m
%
% Reading stuff on race and IQ by Steve Gould, Murray & Herrnstein
%
% Given the population difference in IQ between two populations, what is
% the positive predictive value for an individual?
%
% RTB wrote it, 01 Feb. 2021, NEH Covid-19 quarantine with MNS

% IQ is normally distributed with mean of 100 and s.d. of 15. The generally
% accepted difference between black and white populations is 15.

nSamp = 1000;
muW = 100; muB = 85;
sdW = 15;

% Draws from our populations:
popW = normrnd(muW,sdW,nSamp,1);
popB = normrnd(muB,sdW,nSamp,1);

% This should be equivalent (and is):
% popW = randn(nSamp,1) + 1;
% popB = randn(nSamp,1);

% plot the distributions:
figure
histogram(popW);
hold on
histogram(popB);

auROC = myroc(popB,popW,'b',1);

% Now we ask the question, given a draw from each population, how often
% would we be correct by assigning the higher draw to the W population?
nCorr = sum(popW > popB);
nCorr = nCorr + 0.5*sum(popW == popB);
PPV = nCorr / nSamp;
