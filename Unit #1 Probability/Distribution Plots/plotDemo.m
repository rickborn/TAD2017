% plotDemo.m
%
% Progressively more information rich methods of showing data

%% Load in four mystery distributions
load mysteryDistributions

% Each distribution is a column vector, concatenate to n x 4:
myData = [D1,D2,D3,D4];
myMeans = mean(myData);
%% Simple bar plot of means
figure, bar(myMeans);
hold on
xlabel('Distribution #');
ylabel('Mean value');
%% Add error bars
nObs = size(myData,1);
stdError = std(myData) ./ sqrt(nObs);
errorbar(myMeans,stdError,'r.');

%% Box plot
figure
boxplot(myData);
xlabel('Distribution #');
ylabel('Value');

%% Distribution plot (a.k.a. "violin plot")
figure
distributionPlot(myData,'histOpt',2); % histOpt=2 works better for uniform distributions than the default
title('Distribution or ''Violin'' Plots');

%% Look at different options for the distribution plot function
distributionPlotDemo