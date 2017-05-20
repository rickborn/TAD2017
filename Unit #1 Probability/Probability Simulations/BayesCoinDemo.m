% BayesCoinDemo.m
%
% Based on http://www.nature.com/nmeth/journal/v12/n5/full/nmeth.3368.html
%
% RTB wrote it 04 May 2017 (badly itching eyes)

% Note: I think this could be Ng-ified for an in-class exercise

%% Frequentist approach: We flip a coin 3 times and get 3 heads
% Is the coin fair, or is something fishy going on?

% The frequentist approach is to ask how unlikely our data is under H0
% What is H0? Ans. pHeads = 0.5

% We would typically ask only about the probability of the data under H0
p3H = 1 - binocdf(2,3,0.5); % 0.125 for 1-tailed (Heads only)
% or could use: binopdf(3,3,0.5)

% We might also calculate a 95% confidence interval
[pHat,pCI] = binofit(3,3,0.05);
% pHat = 1; pCI = 0.2924    1.0000

%% Bayesian approach: We flip a coin 3 times and get 3 heads
% Here we consider an entire range of possible hypotheses concerning the
% probability of heads:
% range of possible values for P(heads)
pi0 = 0:0.01:1;
nBins = length(pi0);

% If we have no clue as to what the true value of p(heads) is, we could
% consider all the values to be equally likely: flat prior
% NOTE that this is sort of absurd, since it is physically impossible to
% have a coin that always lands heads unless it has heads on both sides, in
% which case it would be immediately obvious from inspection. Perhaps this
% is a teachable moment about priors.
priorFlat = ones(size(pi0)) / nBins;

likelihood3H = binopdf(3,3,pi0);    % general formulation, could just use pi0.^3
likelihood3H = likelihood3H ./ trapz(likelihood3H); % normalize by area

% calculate posterior and normalize
posterior3H = likelihood3H .* priorFlat;
posterior3H = posterior3H ./ trapz(posterior3H);

% Plot 'em
figure, 
subplot(2,2,1)
plot(pi0,priorFlat,'k-',pi0,likelihood3H,'r-',pi0,posterior3H,'b--');
legend('Prior','Likelihood','Posterior', 'Location','northwest');
xlabel('\pi'); ylabel('P(\pi)');
title('Bayes: flat prior');

%% Calculate some values of interest
% Note: This is an opportunity to point out the richness afforded by having
% the entire posterior distribution as opposed to just a point estimate or
% an interval estimate.

% Expectation:
muPosterior3H = sum(posterior3H .* pi0);    % E(pi)

% The median divides the probability mass in half, so we want to find the
% value of pi0, pi0median, for which trapz(posterior3H <= pi0median) = 0.5
% This is a brute force approximation.
idx = 1;
while trapz(posterior3H(pi0 <= pi0(idx))) <= 0.5
    idx = idx+1;
end
pi0median = pi0(idx-1);

% The most likely value of pi0 is the mode:
piMax = pi0(posterior3H == max(posterior3H));

% What is the probability that the coin is biased towards heads?
% Note that this question doesn't exist for the frequentist
pBiased = trapz(posterior3H(pi0 > 0.5));  
% What are the odds the coin is biased towards heads?
oddsBiased2Heads = pBiased / (1-pBiased);

%% Non-flat priors
% Now, what if we base our prior on having observed 3H and 1T:
prior3H1T = (pi0.^3).*(1-pi0);
prior3H1T = prior3H1T ./ trapz(prior3H1T);

% Then we observe an outcome of 4H:
likelihood4H = binopdf(4,4,pi0);
likelihood4H = likelihood4H ./ trapz(likelihood4H);

posterior3H1T = likelihood4H .* prior3H1T;
posterior3H1T = posterior3H1T ./ trapz(posterior3H1T);

subplot(2,2,2);
plot(pi0,prior3H1T,'k-',pi0,likelihood4H,'r--',pi0,posterior3H1T,'b-');
legend('Prior','Likelihood','Posterior', 'Location','northwest');
xlabel('\pi'); ylabel('P(\pi)');
title('Prior based on 3H1T');

%% Gaussian prior centered around 0.5, observe 3H1T
priorGaussFair = normpdf(pi0,0.5,0.25);
priorGaussFair = priorGaussFair ./ trapz(priorGaussFair);

% Then we observe an outcome of 3H1T:
likelihood3H1T = binopdf(3,4,pi0);
likelihood3H1T = likelihood3H1T ./ trapz(likelihood3H1T);

posterior3H1T = likelihood3H1T .* priorGaussFair;
posterior3H1T = posterior3H1T ./ trapz(posterior3H1T);

subplot(2,2,3);
plot(pi0,priorGaussFair,'k-',pi0,likelihood3H1T,'r-',pi0,posterior3H1T,'b-');
legend('Prior','Likelihood','Posterior', 'Location','northwest');
xlabel('\pi'); ylabel('P(\pi)');
title('Gaussian prior around 0.5');
%% Gaussian prior tightly centered around 0.5, observe 3H1T
priorGaussFair = normpdf(pi0,0.5,0.1);
priorGaussFair = priorGaussFair ./ trapz(priorGaussFair);

% Then we observe an outcome of 3H1T:
likelihood3H1T = binopdf(3,4,pi0);
likelihood3H1T = likelihood3H1T ./ trapz(likelihood3H1T);

posterior3H1T = likelihood3H1T .* priorGaussFair;
posterior3H1T = posterior3H1T ./ trapz(posterior3H1T);

subplot(2,2,4);
plot(pi0,priorGaussFair,'k-',pi0,likelihood3H1T,'r-',pi0,posterior3H1T,'b-');
legend('Prior','Likelihood','Posterior', 'Location','northwest');
xlabel('\pi'); ylabel('P(\pi)');
title('Prior tight near 0.5');