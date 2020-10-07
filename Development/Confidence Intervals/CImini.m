% CImini.m: a quick lesson about calculating confidence intervals
%
% RTB wrote it, 03 October 2020, Just returned to JP from NEH, Trump
% Covid-19+

alpha = 0.05;
n = 100;

% Generate random draws from the standard normal distribution:
a = randn(n,1);

% Use 'normfit' to calculate the mean and 95% CI
[muHat,~,muCI,~] = normfit(a);

% Now let's do it 'by hand':
aMu = mean(a);
aSEM = std(a) / sqrt(n);

% We normally think that our 95% CI would be +/- 1.96 SEMs:
aCIupperNl = aMu + aSEM * norminv(1 - alpha/2);
aCIlowerNl = aMu - aSEM * norminv(1 - alpha/2);

% But when we compare these with the values from 'normfit' we see that they
% are close, but not exact. What gives? Recall that we can use the
% z-distribution when we know the population mean and std dev. But, in this
% case, both are based on the sample. The "price" we pay for this is that
% we have to use the t-distribution, which has slightly fatter tails.
aCIupperT = aMu + aSEM * tinv(1 - alpha/2,n-1);
aCIlowerT = aMu - aSEM * tinv(1 - alpha/2,n-1);

% Now we see that these match exactly with the values from 'normfit'