% interactionPowerDemo.m
%
% RTB wrote it, 22 March 2017

% This is from Andrew Gelman's blog, with translations into MATLAB by RTB

% You need 16 times the sample size to estimate an interaction than to
% estimate a main effect
%
% Posted by Andrew on 15 March 2018, 9:11 am

% Yesterday I shared the following exam question:
% 
% In causal inference, it is often important to study varying treatment
% effects: for example, a treatment could be more effective for men than
% for women, or for healthy than for unhealthy patients. Suppose a study is
% designed to have 80% power to detect a main effect at a 95% confidence
% level. Further suppose that interactions of interest are half the size of
% main effects. What is its power for detecting an interaction, comparing
% men to women (say) in a study that is half men and half women? Suppose
% 1000 studies of this size are performed. How many of the studies would
% you expect to report a statistically significant interaction? Of these,
% what is the expectation of the ratio of estimated effect size to actual
% effect size?
% 
% Here’s the solution:
% 
% If you have 80% power, then the underlying effect size for the main
% effect is 2.8 standard errors from zero. That is, the z-score has a mean
% of 2.8 and standard deviation of 1, and there’s an 80% chance that the
% z-score exceeds 1.96 (in R, pnorm(2.8, 1.96, 1) = 0.8).

% RTB: We can use 'norminv' to calculate the effect size directly:
effectSize = norminv(0.8,1.96,1);   % Ans: 2.8016

% NOTE: 'pnorm' in 'R' is equivalent to 'normcdf' in MATLAB:
myPower = normcdf(2.8,1.96,1);      % Ans: 0.7995


% Now to the interaction. The standard of an interaction is roughly twice
% the standard error of the main effect, as we can see from some simple
% algebra: – The estimate of the main effect is ybar_1 – ybar_2, which has
% standard error sqrt(sigma^2/(N/2) + sigma^2/(N/2)) = 2*sigma/sqrt(N); for
% simplicity I’m assuming a constant variance within groups, which will
% typically be a good approximation for binary data, for example. – The
% estimate of the interaction is (ybar_1 – ybar_2) – (ybar_3 – ybar_4),
% which has standard error sqrt(sigma^2/(N/4) + sigma^2/(N/4) +
% sigma^2/(N/4) + sigma^2/(N/4)) = 4*sigma/sqrt(N). [algebra fixed]
% 
% And, from the statement of the problem, we’ve assumed the interaction is
% half the size of the main effect. So if the main effect is 2.8 on some
% scale with a se of 1, then the interaction is 1.4 with an se of 2, thus
% the z-score of the interaction has a mean of 0.7 and a sd of 1, and the
% probability of seeing a statistically significant effect difference is
% pnorm(0.7, 1.96, 1) = 0.10. That’s right: if you have 80% power to
% estimate the main effect, you have 10% power to estimate the interaction.
probSig = normcdf(0.7,1.96,1);      % Ans: 0.1038

% And 10% power is really bad. It’s worse than it looks. 10% power kinda
% looks like it might be OK; after all, it still represents a 10% chance of
% a win. But that’s not right at all: if you do get “statistical
% significance” in that case, your estimate is a huge overestimate:
% 
% > raw <- rnorm(1e6, .7, 1) 
% > significant <- raw > 1.96 
% > mean(raw[significant]) 
% [1] 2.4 

rawData = normrnd(0.7,1,100000,1);
meanEffectSize = mean(rawData(rawData > 1.96));
% Ans: 2.4349

% So, the 10% of results which do appear to be statistically significant
% give an estimate of 2.4, on average, which is over 3 times higher than
% the true effect.
% 
% Take-home point
% 
% The most important point here, though, has nothing to do with statistical
% significance. It’s just this: Based on some reasonable assumptions
% regarding main effects and interactions, you need 16 times the sample
% size to estimate an interaction than to estimate a main effect.
% 
% And this implies a major, major problem with the usual plan of designing
% a study with a focus on the main effect, maybe even preregistering, and
% then looking to see what shows up in the interactions. Or, even worse,
% designing a study, not finding the anticipated main effect, and then
% using the interactions to bail you out. The problem is not just that this
% sort of analysis is “exploratory”; it’s that these data are a lot noisier
% than you realize, so what you think of as interesting exploratory
% findings could be just a bunch of noise.
% 
% I don’t know if all this in the textbooks, but it should be.