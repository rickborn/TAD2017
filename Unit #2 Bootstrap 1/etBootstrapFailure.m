% etBootstrapFailure.m
%
% Example of p. 81 of E&T where nonparametric bootstrap fails
%
% RTB wrote it, 31 Jan. 2017; snow starting to fall

% Suppose we have data, x1,x2 . . . xn, drawn from a uniform distribution
% on (0,theta). The maximum likelihood estimate, thetaHat, is max(x_n)

% Generate a sample of 50 uniform numers over (0,1)
X = rand(50,1);
n = length(X);
thetaHat = max(X);

% Traditional, non-parametric bootstrap
nBoot = 10000;
maxVals = zeros(nBoot,1);

for k=1:nBoot
    maxVals(k) = max(X(unidrnd(n,n,1)));
end
figure, hist(maxVals,20);
xlabel('Nonparametric'); ylabel('#');
title('E&T fig. 7.11, p. 83, left panel');

% Parametric bootstrap: sample from uniform distribution (0,thetaHat)
maxVals = zeros(nBoot,1);

for k=1:nBoot
    maxVals(k) = max(rand(n,1) .* thetaHat);
end
figure, hist(maxVals,20);
xlabel('Parametric'); ylabel('#');
title('E&T fig. 7.11, p. 83, right panel');

% E&T p. 81: "What goes wrong with the nonparametric bootstrap? The
% difficulty occurs because the empirical distribution function, F_hat, is
% not a good estimate of the true distribution, F, in the extreme tail. . .
% . The nonparametric bootstrap can fail in other examples in which theta
% depends on the smootheness of F."