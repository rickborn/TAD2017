% etBootstrapFailure.m
%
% Example of p. 81 of E&T where nonparametric bootstrap fails
%
% RTB wrote it, 31 Jan. 2017; snow starting to fall

% Suppose we have data, x1,x2 . . . xn, drawn from a uniform distribution
% on (0,theta). The maximum likelihood estimate, thetaHat, is max(x_n)

% Generate a sample of 50 uniform numbers over (0,1)
X = rand(50,1);
n = length(X);
thetaHat = max(X);

% Traditional, non-parametric bootstrap
nBoot = 10000;
maxVals = zeros(nBoot,1);

for k=1:nBoot
    maxVals(k) = max(X(unidrnd(n,n,1)));
end
figure('Name','Bootstrap Failure');
subplot(2,1,1);
histogram(maxVals,20);
hold on
xlabel('Nonparametric'); ylabel('#');
title('E&T fig. 7.11, p. 83, left panel');
ax = axis;
line([thetaHat,thetaHat],[ax(3),ax(4)],'Color','r','LineWidth',2);

% Parametric bootstrap: sample from uniform distribution (0,thetaHat)
maxVals = zeros(nBoot,1);

for k=1:nBoot
    maxVals(k) = max(rand(n,1) .* thetaHat);
end
subplot(2,1,2);
histogram(maxVals,20);
hold on
xlabel('Parametric'); ylabel('#');
title('E&T fig. 7.11, p. 83, right panel');
ax = axis;
line([thetaHat,thetaHat],[ax(3),ax(4)],'Color','r','LineWidth',2);

% E&T p. 81: "What goes wrong with the nonparametric bootstrap? The
% difficulty occurs because the empirical distribution function, F_hat, is
% not a good estimate of the true distribution, F, in the extreme tail. . .
% . The nonparametric bootstrap can fail in other examples in which theta
% depends on the smootheness of F."

% From the top histogram, we see that our nonparametric bootstrap replicate
% was exactly equal to thetaHat 62% of the time. This is just the
% probability of any given value in our original sample of being included
% in any given bootstrap sample: 1 - (1 - 1/n)^n, which approaches 1 - 1/e
% as n approaches infinity, which is 1 - exp(-1) = 0.6321.