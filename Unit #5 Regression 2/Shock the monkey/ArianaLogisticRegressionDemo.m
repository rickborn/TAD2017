% ArianaLogisticRegressionDemo.m
%
% For stimulus detection paradigm to be used in the micro-coil experiments
% Plot detection data and fit a logistic regression model.
%
% RTB wrote it, 25 April 2018

%% Read in the data 
cd 'Z:\Born Lab\Ariana\starborn'
startup
s = detectg('/');

% Select the file you want to analyze. (for example the str file in
% Ariana’s folder on research.files, in particular:
% data/2018/04/24/ariana_11.18.46_gaussdetection.str.

D = gather_responses(s);

% D is a matrix with one row for each trial. The columns have the different
% data values in them:
% col 1: stimulus detected (1) or not detected (0), (-1 for early response)
% col 2: reaction time (ms)
% col 3: gain of stimulus (measure of luminance)
% cols 4 and 5: x and y position of stimulus (°).

% Manifest constants for column labels
DETECT = 1;
RT = 2;
LUM = 3;
STIMX = 4;
STIMY = 5;
STIMECC = 6;

%% Calculate a column for eccentricity

% remove 'early' responses
selL = D(:,DETECT) == -1;
D(selL,:) = [];

% Calculate stimulus eccentricity from X and Y locations
D(:,STIMECC) = sqrt(D(:,STIMX).^2 + D(:,STIMY).^2);

%% Plot the data
% We want a plot like those in fig. 4 of Salzman et al. 1992:
% x-axis is the signed correlation value; y-axis is proportion of preferred
% decisions.

% We first need to find all stimulus conditions:
allDataProp = sortrows(unique([D(:,LUM),D(:,STIMECC)],'rows'),[1,2]);
nConds = length(allDataProp);
allLums = unique(D(:,LUM),'rows');
allEccs = unique(D(:,STIMECC),'rows');

% Add a column of zeros to allDataProp to hold the calcualted proportion;
% (plus two more for confidence intervals)
allDataProp = [allDataProp,zeros(nConds,3)];
% myErr = 0.32;   % for error bars of the 68% CI (= SEM)
myErr = 0.05;   % for error bars of the 95% CI
for k = 1:nConds
    thisCond = find(D(:,LUM) == allDataProp(k,1) & D(:,STIMECC) == allDataProp(k,2));
    % allDataProp(k,3) = sum(ds.PDchoice(thisCond)) / length(thisCond);
    [pHat,pCI] = binofit(sum(D(thisCond,DETECT)),length(thisCond),myErr);
    allDataProp(k,3) = pHat;
    allDataProp(k,[4,5]) = pCI;
end

% Convert to values appropriate for errorbar function:
allDataProp(:,4) = allDataProp(:,3) - allDataProp(:,4); % lower error bar
allDataProp(:,5) = allDataProp(:,5) - allDataProp(:,3); % upper error bar

figure, hold
for k = 1:length(allEccs)
    thisData = allDataProp(allDataProp(:,2) == allEccs(k),:);
    hP(k) = errorbar(thisData(:,1),thisData(:,3),...
    thisData(:,4),thisData(:,5),'*-');
end
% Now plot: stim trials in red; no-stim trials in black
xlabel('Luminance'); ylabel('Proportion detected');
set(gca,'FontSize',14); % make text larger
legend(hP,num2str(allEccs,2),'Location','SouthEast');

%% Fit model using glmfit

% ln(P/1-P) = b0 + b1*lum + b2*ecc
[b,dev,stats] = glmfit([D(:,LUM),D(:,STIMECC)],D(:,DETECT),'binomial','link','logit');

T = table(b,stats.se,stats.p,'VariableNames',{'Beta','SE','p'});
disp(T);

%% Plot the model fits

lum = min(D(:,LUM)):0.01:max(D(:,LUM));
const = ones(size(lum));
pDetect = 1 ./ (1 + exp(-(b(1).*const + b(2).*lum)));
plot(lum,pDetect,'k-','LineWidth',2);
