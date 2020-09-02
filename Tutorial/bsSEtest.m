% bsSEtest.m

% shuffle the random number generator
rng shuffle

% 100 draws from std nl distribution
A = normrnd(0,1,100,1);
figure, histogram(A)

[muhat,stdhat,muci,stdci] = normfit(A);
sehat = stdhat / sqrt(length(A));

% estimate the standard error by bootstrapping
nBoot = 100000;
allMeans = zeros(nBoot,1);
for k = 1:nBoot
    % re-sample from A with replacement:
    allMeans(k) = mean(A(unidrnd(length(A),length(A),1)));
end
figure, histogram(allMeans);