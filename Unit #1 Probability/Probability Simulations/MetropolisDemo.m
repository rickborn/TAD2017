function [allTh, allProp] = MetropolisDemo(nSims)

% MetropolisDemo.m: simple algorithm for Markov Chain Monte Carlo (MCMC)
%
% A 1D instantiation of the Metropolis Algorithm
% from Chapter 7 of Kruschke, beginning on p. 146 (7.2.1)
%
% Metropolis,N., Rosenbluth,A.W., Rosenbluth,M.N., Teller,A.H.,&
% Teller,E.(1953). Equations of state calculations by fast computing
% machines. Journal of Chemical Physics, 21, 1087-1091.

% Note on terminology. Any process that samples random values from some
% distribution is called a Monte Carlo process (von Neumann coined the
% term). Any random walk in which each step is independent of the steps
% before the current position (i.e. no memory) is called a (1st order)
% Markov process, and the succession of such steps is called a Markov
% chain. So the Metropolis algorithm and Gibbs sampling are just particular
% examples of a Markov Chain Monte Carlo (MCMC) process.

% The analogy: "Suppose an elected politician lives on a long chain of
% islands. He is constantly traveling from island to island, wanting to
% stay in the public eye. At the end of a grueling day of photo
% opportunities and fundraising, he has to decided whether to (i) stay on
% the current island, (ii) move to the adjacent island to the west, or
% (iii) move to the adjacent island to the east. His goal is to visit all
% the islands proportionally to their relative population, so that he
% spends the most time on the most populated islands, and proportionally
% less time on the less populated islands."

% Heuristic of the politician for visiting each island. Flip a fair coin to
% determine whether to move to the east (right) or west (left). If the
% proposed island has a larger population, then he definitely goes to the
% posposed island. If not, then he goes to the proposed island with
% probability of P_proposed / P_current.

% The plot re-capitulates Figure 7.2 of Kruschke

if nargin < 1
    nSims = 10000;
end
rng default

th = 0:8;   % which island
pTh = [0,round(rand(1,7).*7),0];    % population of each island

% The original algorithm fails with something like this, because we are
% only able to jump 1 island, which means having a 0 in the middle acts as
% an insurmountable barrier. It can be fixed by allowing for larger jumps.
%pTh = [4,3,2,0,2,3,4];

figure;
subplot(3,1,1);
bar(th(2:8),pTh(2:8));
xlabel('\theta');
ylabel('P(\theta)');
title('Target distribution');

allTh = zeros(nSims,1);
allProp = zeros(nSims,1);
thisTh = 4;     % start on island #4; could be random
allTh(1) = thisTh; 
allProp(1) = thisTh;

for k = 2:nSims
    % We propose to jump to a random island: "proposal distribution"
    propTh = thisTh + randsample([-2,-1,1,2],1);
%     while propTh < 1 || propTh > 7
%         propTh = thisTh + randsample([-3,-2,-1,1,2,3],1);
%     end
    allProp(k) = propTh;
    % If the population is greater on the proposed island, we go there:
    if pTh(th==propTh) > pTh(th==thisTh)
        allTh(k) = propTh;
    % If not, we go to the proposed island with probability:
    % population_proposed / population_current
    else
        if rand < (pTh(th==propTh)/pTh(th==thisTh))
            allTh(k) = propTh;
        else
            allTh(k) = thisTh;
        end
    end
    thisTh = allTh(k);
end

if nSims < 500
    n = nSims;
else
    n=500;
end
subplot(3,1,2);
semilogy(allTh(1:n),1:n,'bo-');
xlabel('\theta');
ylabel('Time step');
title('Biased Random walk');

subplot(3,1,3);
xBins = 1:7;
allCounts = hist(allTh,xBins);
bar(xBins,allCounts);
ax = axis;
axis([ax(1),ax(2),ax(3),max(allCounts)]);
xlabel('\theta');
ylabel('Frequency');
title(['Sampled distribution, n = ' num2str(nSims)]);

% Post hoc thoughts. To make this work, we only need to be able to do 3
% things:
%
% 1) generate a random value from the "proposal distribution"
% 2) evaluate the target distribution at any proposed position
% 3) generate a random value from a uniform distribution