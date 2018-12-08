% suspiciousDataDemo.m
%
% RTB wrote it, 05 Dec. 2018, post-triage reward
%
% inspired by http://datacolada.org/74#identifier_4_2860
% data from: https://osf.io/k8tv2/

% cd 'C:\usr\rick\doc\Committees\PIN\PIN Director\Courses\Stats\TAD\TAD Code\Development\Hand Sanitizer Study'

% Data are daily hand sanitizer usage for two groups over 40 days. The
% first 20 days were to establish baseline; the treatment starts on day 21
%
% In Experiment 2, sanitizer use was measured daily for 40 participants for
% 40 days (20 days of baseline, 20 of treatment), all with a scale
% sensitive to 100th of a gram. Recall that the manipulation was done at
% the room level. This figure, which was in the original article, shows the
% daily average use of sanitizer across the two rooms.

%% Read in data and plot it
[num,txt,raw] = xlsread('Data from all studies.xlsx',4);

% vector for selecting experimental vs. control groups
expID = logical(num(:,2));
% summary stats for the two groups
muExp = mean(num(expID,3:end));
semExp = std(num(expID,3:end)) / sqrt(sum(~expID));
muCtrl = mean(num(~expID,3:end));
semCtrl = std(num(~expID,3:end)) / sqrt(sum(~expID));

% plot the two groups
xDays = 1:40;
figure('position',[50 50 1200 600]);
errorbar(xDays,muExp,semExp,'o-','LineWidth',2);
hold on
errorbar(xDays,muCtrl,semCtrl,'o-','LineWidth',2);
xlabel('Day of Experiment');
ylabel('Mean grams of hand sanitizer used');
title('Li, Sun & Chen: Decoy effect as a nudge');

% treatment starts on day 21; draw a line
ax = axis;
h1 = line([20.5,20.5],[ax(3),ax(4)],'Color','k','LineStyle','--');

legend('Spray bottle + soak basin (decoy)','Spray bottle only (no decoy)',...
    'Start of Rx','Location','northwest');
set(gca,'FontSize',14);

%% Implausibly similar means

% Anything suspicious? During the "before treatment" period, the two groups
% are surprisingly similar. How correlated are the means?
rhoExp = corrcoef(muExp(1:20),muCtrl(1:20));
rhoExp = rhoExp(2);
%display(rhoExp);

muSqDiffExp = mean((muExp(1:20)- muCtrl(1:20)).^2);

%% Permutation test

% How surprising is this? Let's do a permutation test on the pre-rx data:
allData = num(:,3:22);

% We want to randomly assign each subject to the rx vs. ctrl group
nPerm = 100000;
allRho = zeros(nPerm,1);
allDiff = zeros(nPerm,1);
for k = 1:nPerm
    shuffledData = allData(randperm(40),:);
    muShuffled1 = mean(shuffledData(1:20,:));
    muShuffled2 = mean(shuffledData(21:end,:));
    r = corrcoef(muShuffled1,muShuffled2);
    allRho(k) = r(2);
    allDiff(k) = mean((muShuffled1 - muShuffled2).^2);
end

figure('position',[50 50 1200 600]);

subplot(1,2,1);
histogram(allDiff);
hold on;
xlabel('Mean squared difference of means');
ylabel('Frequency');
title(['Study 2: Similarity of means across rooms in ' num2str(nPerm) ' permuted samples']);
ax = axis;
line([muSqDiffExp,muSqDiffExp],[ax(3),ax(4)],'Color','r');

pVal1 = sum(allDiff <= muSqDiffExp) / nPerm;
text(40,(ax(4)-ax(3))/1.5,['p-value: ' num2str(pVal1)]);

subplot(1,2,2);
histogram(allRho);
hold on
xlabel('Correlation of means across days');
ylabel('Frequency');
title(['Study 2: Correlation of means across rooms in ' num2str(nPerm) ' permuted samples']);
ax = axis;
line([rhoExp,rhoExp],[ax(3),ax(4)],'Color','r');

pVal2 = sum(allRho >= rhoExp) / nPerm;
% display(pVal);
text(0.3,(ax(4)-ax(3))/1.5,['p-value: ' num2str(pVal2)]);

%% Violation of Benford's Law?

% see: https://en.wikipedia.org/wiki/Benford%27s_law

% Make a histogram of the frequency of each of the last digits. All values
% are recorded to the nearest 100th of a gram, so this is easy:
allNum = num(:,3:end);
lastDigits = mod(round(allNum(:) .* 100), 10);
figure('position',[50 50 1200 600]);
% observed counts per bin:
subplot(1,2,1);
O = hist(lastDigits,[0:9]);
bar([0:9],O);
xlabel('Last digit (1/100th gram)');
ylabel('Frequency');
title('A violation of Benford''s Law?');

% prediction for uniformity:
nPerBinUni = numel(allNum) / 10;
ax = axis;
axis([-1,10,ax(3),ax(4)]);
h = line([-1,10],[nPerBinUni,nPerBinUni],'Color','r');

% Chi-squared test for uniformity.
E = repmat(nPerBinUni,size(O));
chi2 = sum(((E-O).^2) ./ E);
df = length(E) - 1;
pChi2 = 1 - chi2cdf(chi2,df);

tStr = sprintf('\\chi^2(%d)=%.2f, p=%.3e',df,chi2,pChi2);
text(1,200,tStr);
%text(1,200,['p-value: ' num2str(pChi2)]);

%% Intuition about H0 distribution of chi2 statistic

% Generate the null distribution under expectation of uniformity
nSim = 100000;
nBins = length(E);
allChiSq = zeros(nSim,1);
N = numel(allNum);
% The expected distribution under H0 (i.e. exactly 160 in each bin)
E = repmat(nPerBinUni,size(O));

for k = 1:nSim
    simData = unidrnd(nBins,N,1) - 1;
    O = hist(simData,[0:9]);
    allChiSq(k) = sum(((E-O).^2) ./ E);
end

subplot(1,2,2);
% histogram(allChiSq);
yyaxis left
[nCts,xBins] = hist(allChiSq,35);
bar(xBins,nCts);
xlabel('\chi^2 Statistic');
ylabel('Frequency');
title('Simulating H0 of the relevant \chi^2 statistic');

% ToDo: yyplot and make a curve of the chi2cdf as well
yyaxis right
yProb = chi2pdf(xBins,df);
plot(xBins,yProb,'.-','Color',[0.8,0.3,0],'LineWidth',2);
ylabel('\chi^2 Probability');

% Draw a line for the value of chisq we actually obtained
% ax = axis;
% line([chi2,chi2],[ax(3),ax(4)],'Color','k');
