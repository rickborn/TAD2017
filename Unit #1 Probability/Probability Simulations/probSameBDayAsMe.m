function [pSim,pExact] = probSameBDayAsMe(n)

% Assume there are 70 students in this class. What is the probability (in
% percent) that one or more other students in this class has the same
% birthday as you? Round to one decimal point. (For this and all other
% parts of this exercise, you can ignore leap years).

% By probability calculus, it is just 1 minus the probability that no other
% student has the same birthday: 1 - (364/365)^n-1 = 0.1725 (for n = 70)
pExact = 1 - (364/365)^(n-1);

% By simulation:
nSims = 1000000;     % # of simulations to run
R = 365;             % # of days in a year (ignore leap years)

% Arbitrarily pick a day for 'my' birthday
myBDay = randi(R);

allMatches = sum(randi(R,n-1,nSims) == myBDay);
pSim = sum(allMatches > 0) / nSims;