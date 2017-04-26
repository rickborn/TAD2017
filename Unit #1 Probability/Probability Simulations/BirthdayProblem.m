function [allP] = BirthdayProblem(K,allN)
%
% Given a group of N people, what is the probability that any K of them
% have the same birthday?
%
% Inputs:
% - K: number with same birthday (default = 2)
% - allN: a range of crowd sizes to test (default = 1:100)
%
% Outputs:
% - allP: probability corresponding to each crowd size
%
% See also: ProbHist.m
% See also: http://en.wikipedia.org/wiki/Birthday_problem
%
% RTB wrote it, 09 August 2013

if nargin < 2, allN = 1:100; end
if nargin < 1, K = 2; end

R = 365; nSims = 10000;

allP = zeros(size(allN));

for j = 1:length(allN)
    allP(j) = sum(any(hist(randi(R,allN(j),nSims),1:R) >= K)) / nSims;
end

figure
plot(allN,allP,'r-');
grid on
xlabel('# of people in group');
yStr = sprintf('probability that %d or more share a birthday',K);
ylabel(yStr);
title('The Birthday Problem');