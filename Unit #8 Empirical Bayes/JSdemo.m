% JSdemo.m: Demo of the James-Stein estimator with batting averages
%
% RTB wrote it 05 Oct. 2016, procrastinating from reading T32s

%% Read in batting average data
ds = dataset('xlsfile','BattingAverages.xlsx');

%% Calculate and plot the maximum likelihood estimator

% We first need to convert the text batting averages ('18/45') to numbers
muMLE = zeros(length(ds),1);
for k = 1:length(ds)
    muMLE(k) = eval(ds.hitAB{k});
end

x = [1:length(ds)];     % x values for plots
p1 = plot(x,muMLE,'k*');
hold on
set(gca,'XTick',[1:18],'XTickLabel',ds.Player,'XTickLabelRotation',60);
ylabel('Batting Average');

% line for the grand mean across players
hl=line([1,length(ds)],[mean(muMLE),mean(muMLE)]);
set(hl,'Color','k','LineStyle','--');

%% Calculate and plot the James-Stein estimator
muJS = JSestimator(muMLE,45);
p2 = plot(x,muJS,'ro');

%% Plot the season batting average--this is the ground truth value
p3 = plot(x,ds.SeasonBA,'bs');
set(p3,'MarkerFaceColor','b');

% add a legend
legend([p1,p2,p3],'mu^M^L^E','mu^J^S','mu');

%% Compute the ratio of prediction errors for the two estimators

RPE = sum((muJS-ds.SeasonBA).^2) / sum((muMLE - ds.SeasonBA).^2);
tStr = sprintf('Ratio of prediction errors: mu^J^S/mu^M^L^E = %0.2f', RPE);
title(tStr);