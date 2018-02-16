% BayesNailCoinDemo.m
%
% RTB wrote it, 10 January 2018; end of bitter cold snap

% This exercise is derived from Chapter 11 (pp. 315-17) of John Kruscheke's
% book on Bayesian Data Analyis

% The Bayesian statistician starts theanalysis with anexpression of the
% prior knowledge. We know from prior experience that the narrow-headed
% nail is biased to show tails, so we experess that knowledge in a prior
% distribution. In our present fictional example involving a nail, suppose
% that we represent our prior beliefs by a fictitious previuos sample that
% had 95% tails in a sample size of 20. That translates into a beta(?|2,20)
% prior distribution if the "proto-prior," before the fictional data, was
% beta(?|1,1).
grain = 0.001;
pH = 0:grain:1;  % range of hypotheses for p(Heads)
priorNail = betapdf(pH,2,20);

% On the other hand, if we have prior knowledge that the object is fair,
% such as a coin,then the prior distribution is different. Suppose that we
% represent our prior beliefs by a fictitious previous sample that had 50%
% tails in a sample size of 20. That translates to a beta(?|11,11) prior
% distribution if the "proto-prior," before the fictional data, was
% beta(?|1,1).
priorCoin = betapdf(pH,11,11);

% Now we do an experiment in which we toss of nail/coin 24 times and get 7
% heads. The likelihood is straightforward:
like1 = binopdf(7,24,pH);

% We can calculate our posterior in one of two ways. One is just to
% multiply the prior x likelihood; the other is to know that it is just
% another, updated beta distribution:
% for the nail: beta(?|2+7,20+17)
% for the coin: beta(?|11+7,11+17)
postNail = betapdf(pH,9,37);
postCoin = betapdf(pH,18,28);
% postNail = priorNail .* like1;
% postCoin = priorCoin .* like1;

% plot these
figure

% prior for the nail
subplot(3,2,1)
bar(pH,priorNail);
ax=axis;
axis([0,1,ax(3),ax(4)]);
xlabel('\theta');
ylabel('dbeta(\theta|2,20)')
title('Prior (nail)');

% prior for the coin
subplot(3,2,2)
bar(pH,priorCoin);
ax=axis;
axis([0,1,ax(3),ax(4)]);
xlabel('\theta');
ylabel('dbeta(\theta|11,11)')
title('Prior (coin)');

% likelihood for the nail
subplot(3,2,3)
bar(pH,like1);
ax=axis;
axis([0,1,ax(3),ax(4)]);
xlabel('\theta');
ylabel('p(D|\theta)')
title('Likelihood (nail)');

% likelihood for the coin
subplot(3,2,4)
bar(pH,like1);
ax=axis;
axis([0,1,ax(3),ax(4)]);
xlabel('\theta');
ylabel('p(D|\theta)')
title('Likelihood (coin)');

% posterior for the nail
subplot(3,2,5)
bar(pH,postNail);
ax=axis;
axis([0,1,ax(3),ax(4)]);
xlabel('\theta');
ylabel('dbeta(\theta|9,37)')
title('Posterior (nail)');

% posterior for the coin
subplot(3,2,6)
bar(pH,postCoin);
ax=axis;
axis([0,1,ax(3),ax(4)]);
xlabel('\theta');
ylabel('dbeta(\theta|18,28)')
title('Posterior (coin)');

% Find the 95% highest density interval for each posterior distribution:
% First, find the maximum of each distribution so that we can move outwards
% symmetrically:
maxPnail = pH(postNail == max(postNail));
maxPcoin = pH(postCoin == max(postCoin));

% But we actually want the index, not the value:
maxPnailNdx = find(postNail == max(postNail));
maxPcoinNdx = find(postCoin == max(postCoin));

% Now calculate the critical value for 95% of the area:
myAlpha = 0.05;
critValNail = (1 - myAlpha) * trapz(postNail);
critValCoin = (1 - myAlpha) * trapz(postCoin);

% Now we just start in the middle and work our way out symmetrically until
% we surpass 95% of the area:
myVal = 0;
ctr = 0;
while myVal < critValNail
    ctr = ctr+1;
    myVal = trapz(postNail(maxPnailNdx-ctr:maxPnailNdx+ctr));
end
% HDI = "highest density interval"
HDI95nail = [pH(maxPnailNdx - ctr),pH(maxPnailNdx + ctr)];

myVal = 0;
ctr = 0;
while myVal < critValCoin
    ctr = ctr+1;
    myVal = trapz(postCoin(maxPcoinNdx-ctr:maxPcoinNdx+ctr));
end
HDI95coin = [pH(maxPcoinNdx - ctr),pH(maxPcoinNdx + ctr)];

% NOTE: This method is not perfect, since it moves out symmetrically in the
% index values, NOT in area. Insofar as the posterior distributions are not
% symmetrical, the HDIs will be shifted a bit.
