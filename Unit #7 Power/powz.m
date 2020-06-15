% powz.m: power based on a z-score example
%
% motivated by Tversky & Kahneman, "Belief in the law of small numbers,"
% Psychol. Bull. 1971

% "Suppose you have run an experiment on 20 subjects, and have obtained a
% significant result which confirms your theory (z = 2.23, p < .05,
% two-tailed). You now have cause to run an additional group of 10
% subjects. What do you think the probability is that the results will be
% significant, by a one-tailed test, separately for this group?"

% The critical value of z for a one-tailed test is 1.645:
% norminv(1-0.05) = 1.6449

nSim = 100000;
zReal = 2.23;
nOri = 20;
nRep = 10;
allZ = zeros(nSim,1);
zCrit = norminv(1-0.05);

for k = 1:nSim
    % I think we need to adjust sigma for the different n's
    x = normrnd(zReal,sqrt(nOri/nRep),[nRep,1]);
    allZ(k) = mean(x) / std(x);
end

p = sum(allZ >= zCrit) / nSim;

% I get a probability of just under 0.5