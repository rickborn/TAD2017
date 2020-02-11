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
semExp = std(num(expID,3:end)) / sqrt(sum(expID));
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
% We shuffle *rows* to preserve each S's values across days
nPerm = 10000;
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
nSim = 10000;
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

%% Number bunching (see also my 'bunchtest.m')

% see: http://datacolada.org/77#identifier_0_4081

% In this post I show how one can analyze the frequency with which values
% get repeated within a dataset – what I call "number-bunching" – to
% statistically identify whether the data were likely tampered with. Unlike
% Benford's law, and its generalizations, this approach examines the
% entire number at once, not only the first or last digit.

% The formula for average-frequency, AF, is AF=sum(fi^2)/sum(fi) where fi is
% the frequency of each distinct value i in the data. This metric is
% similar in spirit and formulation to entropyThe formula for
% average-frequency, AF, is AF=sum(fi2)/sum(fi) where fi is the frequency
% of each distinct value i in the data. This metric is similar in spirit
% and formulation to entropy

% Hand sanitizer usage for study #2:
[num,~,~] = xlsread('Data from all studies.xlsx',4);

% Hand sanitizer usage for study #3:
%[num,~,~] = xlsread('Data from all studies.xlsx',7);

allVals = num(:,3:end);      % all of the values in experiment #2
allVals = allVals(:);

% remove missing values (NaN's)
if sum(isnan(allVals))
    allVals(isnan(allVals)) = [];
end

% Find the frequency of each unique value:
uVals = unique(allVals);

% Count up the number of times each distinct value is repeated
fVals = zeros(size(uVals));
for k = 1:length(uVals)
    fVals(k) = sum(allVals == uVals(k));
end
AF = sum(fVals.^2) / sum(fVals);

% Expected bunching:
% Because sanitizer use was measured to 1/100th of a gram, we can generate
% expected levels of number bunching by splitting gram values into an
% integer and decimal part, and then shuffling the decimal parts across
% observations. This effectively treats the decimal portion as random. This
% preserves the overall distribution of weights, as well as the exact
% frequency of every observed decimal value.  All we are assuming here is
% that the (recorded) weight in grams is independent of its decimal portion
% (any observed number of grams could be associated with any observed
% fraction of grams). This test is thus independent of the already shown
% evidence that these data fail the uniform last digit test, for the
% (non-uniform) frequency of last-digits is maintained in the bootstrapped
% samples.

intVals = floor(allVals);
decVals = allVals - intVals;
decVals = round(decVals .* 100) ./ 100;
nVals = numel(allVals);

nSims = 10000;
allAF = zeros(nSims,1);
for j = 1:nSims
    % shuffle the decimal parts and add back to the integer parts
    shuffledDecimals = decVals(randperm(nVals));
    newVals = shuffledDecimals + intVals;
    
    uVals = unique(newVals);
    fVals = zeros(size(uVals));
    for k = 1:length(uVals)
        fVals(k) = sum(newVals == uVals(k));
    end
    allAF(j) = sum(fVals.^2) / sum(fVals);
end
figure, histogram(allAF);
xlabel('Average frequency of observed values');
ylabel('# of simulations');
ax = axis;
line([AF,AF],[ax(3),ax(4)],'Color','r');

% calculate a p-value
pVal = sum(allAF >= AF) ./ nSims;
if pVal == 0
    pVal = 1 / (nSims+1);
end

%% Placebo tests

% The second "placebo" dataset consists of a montecarlo simulation,
% generating data from the normal distribution matching the N, mean, and SD
% from Study 3. (i.e. We expect this to not exhibit unexpected bunching)
n3 = 3234;
mu3 = 54.1;
sd3 = 14.87;
mcVals = normrnd(mu3,sd3,[n3,1]);
mcVals = round(mcVals .* 100) ./ 100;   % round to nearest 100th gram
[p,allSim] = bunchtest(mcVals,2,10000,1);

% For other examples:
% 1) non-suspicious data: BroodCarcassMass.xlsx
%    https://datadryad.org/stash/dataset/doi:10.5061/dryad.3kh41
%
% [num,~,~] = xlsread('BroodCarcassMass.xlsx');
% D = num(:,1);
% [p,allSim] = bunchtest(D,2,10000,1);
% % passes easily
%
% 2) data from a retracted paper by Laskowski et al.
%    https://datadryad.org/stash/dataset/doi:10.5061/dryad.33f0n
% see also: https://laskowskilab.faculty.ucdavis.edu/2020/01/29/retractions/
%
% [num,~,~] = xlsread('Laskowski et al_social niche disruption_DATA.xlsx');
% D = num(:,7:8);
% [p,allSim] = bunchtest(D,2,10000,1);
% % fails badly
%
% 3) data from a paper being investigated for retraction by Grinsted et al.
%   https://datadryad.org/stash/dataset/doi:10.5061/dryad.nd779
%
% [num,~,~] = xlsread('Grinsted et al. 2013 Proceedings B online data.xlsx',6);
% D = num(:,12);
% [p,allSim] = bunchtest(D,2,10000,1);
% % fails badly