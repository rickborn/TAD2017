function [] = CLTdemo(R,nSum,nSim,whichDist)

% CLTdemo.m: Demonstration of the Central Limit Theorem
%
% [] = CLTdemo(R,nSum,nSim,whichDist);
% e.g. CLTdemo(5,[1:5],10000,'unidrnd');
%
% The normal distribution is what you get when you add up a bunch of random
% events.
%
% Inputs:
% - R, the parameter for the discrete uniform distribution (default = 5)
% - nSim, the number of draws to make (default = 10000)
% - nSum, the number of random values to sum (default = [1:5])
% - whichDist, the random distribution to use (default = discrete uniform)
%
% see chapter 7 in Vickers "What is a p-value anyway?"
% RTB wrote it, fall of 2012
% Rebecca's test 5/23/17

% Other distributions to try: exprnd, chi2rnd, betarnd

if nargin < 4, whichDist = 'unidrnd'; end
if nargin < 3, nSim = 10000; end
if nargin < 2, nSum = [1:5]; end
if nargin < 1, R = 5; end

% Start with draws from a non-normal distribution: e.g. uniform
% Show that the sum of several random things approaches normal:
for thisSum = nSum
    fStr = sprintf('allSums = sum(%s(%d,[%d,%d]),1);',whichDist,R,thisSum,nSim);
    eval(fStr);
    %allSums = sum(unidrnd(R,[thisSum nSim]),1);
    xBins = min(allSums):max(allSums);
    myHist = hist(allSums, xBins);
    normHist = myHist ./ sum(myHist);
    figure, bar(xBins, normHist); hold on;
    tStr = sprintf('%s R = %d, nSim = %d',whichDist,R,nSim); title(tStr);
    xStr = sprintf('Sum of %d random draws',thisSum); xlabel(xStr);
    ylabel('Probability');
    
    muSum = mean(allSums);
    sdSum = std(allSums);
    x = min(allSums)-R:0.01:max(allSums)+R;
    y = normpdf(x,muSum,sdSum);
    plot(x,y,'r-');
end

% Another fun variant: sum up a bunch of DIFFERENT distributions
nReps = 10;                      % number of each type to sum
nSim = 10000;
a = round(rand(nReps,nSim));    % 0s & 1s
figure, hist(a',[0,1]);
tempSum = sum(a);

xBins = min(tempSum):max(tempSum);
myHist = hist(tempSum, xBins);
normHist = myHist ./ sum(myHist);
figure, bar(xBins, normHist); hold on;

muSum = mean(tempSum);
sdSum = std(tempSum);
x = min(tempSum):0.01:max(tempSum);
y = normpdf(x,muSum,sdSum);
plot(x,y,'r-');

b = binornd(5,0.4,nReps,nSim);
c = exprnd(2,nReps,nSim);
figure, hist(sum([a;b;c]));