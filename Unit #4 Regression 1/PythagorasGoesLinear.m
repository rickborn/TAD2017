% PythagorasGoesLinear.m
%
% from p. 160-1 of Gelman & Nolan's "Teaching statistics: A bag of tricks"
%
% From p. 161: "What is the point of this example? At one level, it shows
% the power of multiple regression--even when the data come from an
% entirely different model, the regression can fit well. There is also a
% cautionary message, showing the limitations of any purely data-analytic
% method for finding true underlying relations. As we tell the students, if
% Pythagoras knew about multiple regression, he might never have discovered
% his famous theorem!"
%
% IMHO this is also a good lesson on the importance of regression
% diagnostics: the residuals-vs-fitted and qqplots clearly reveal that our
% assumptions for linear regression are not met.
%
% RTB wrote it, 26 Sept. 2017

%% Generate a Pythagorean data set

nSamp = 50;
% x1 = unidrnd(10,nSamp,1);
% x2 = unidrnd(20,nSamp,1);
x1 = rand(nSamp,1) .* 10;
x2 = rand(nSamp,1) .* 20;
y = sqrt(x1.^2 + x2.^2);

%% Model it with linear regression
const = ones(length(y),1);
[betaFit,betaCI,resid,residInt,stats] = regress(y,[const,x1,x2]);

%% Plot residuals vs. fitted

% model predictions:
yPred = betaFit(1) + betaFit(2).*x1 + betaFit(3).*x2;

figure, plot(yPred,resid,'ko');
xlabel('Predicted values');
ylabel('Residuals');
tStr = sprintf('Linear prediction R^2 = %.3f', stats(1));
title(tStr);

%% Q-Q Plot of residuals

figure, qqplot(resid);