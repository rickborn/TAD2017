% plotDemo.m
%
% Progressively more information rich methods of showing data

%% Load in four mystery distributions
load mysteryDistributions

% Each distribution is a column vector, concatenate to n x 4:
myData = [D1,D2,D3,D4];
myMeans = mean(myData);
%% Simple bar plot of means
figure, subplot(1,3,1);
bar(myMeans);
hold on
xlabel('Distribution #');
ylabel('Mean value');
%% Add error bars
nObs = size(myData,1);
stdError = std(myData) ./ sqrt(nObs);
errorbar(myMeans,stdError,'r.');
title('Bar Plots');
%% Box plot
subplot(1,3,2);
boxplot(myData);
xlabel('Distribution #');
ylabel('Value');
title('Box Plots');
%% Distribution plot (a.k.a. "violin plot")
subplot(1,3,3);
distributionPlot(myData,'histOpt',2); % histOpt=2 works better for uniform distributions than the default
title('Violin Plots');
xlabel('Distribution #');
ylabel('Value');
%% Look at different options for the distribution plot function
%distributionPlotDemo