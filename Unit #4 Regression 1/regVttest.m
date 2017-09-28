% regVttest.m
%
% Getting the same answer with regression and a t-test
%
% RTB wrote it, 22 Sept. 2017 (hurricane Jose rainy equinox)

% Using data from one of Brian Healy's practicums: bpf.xls
%
% Each row is one patient
% Each column is a variable:
% col 1: age in years
% col 2: bpf ( = 'brain parenchymal fraction')
% col 3: gender (M,F)
% col 4: male (1 = male, 0 = female)

%% Read in data
ds = dataset('xlsfile','bpf');

male = logical(ds.male);
nMale = sum(male);
nFemale = sum(~male);

%% Make a few plots

figure, gscatter(ds.age,ds.bpf,ds.gender,'br','xo');
xlabel('age (yrs)'); ylabel('brain parenchymal factor');

% figure, scatterhist(ds.age,ds.bpf,'Group',ds.gender,'Location','Northwest',...
%         'Direction','out');
% xlabel('age (yrs)'); ylabel('brain parenchymal factor');

figure, boxplot(ds.bpf,male);
hold on

ylabel('brain parenchymal factor');
set(gca,'XTickLabel',['F';'M'])

%% Superimpose raw data on the box plots

% To do this, we need to know the numeric values associated with the two
% groups. This can be obtained using the 'gca' command
xVals = get(gca,'XTick');

semMale = std(ds.bpf(male)) / sqrt(sum(ds.bpf(male)));
semFemale = std(ds.bpf(~male)) / sqrt(sum(ds.bpf(~male)));
he = errorbar(xVals,[mean(ds.bpf(~male)),mean(ds.bpf(male))],[semFemale,semMale],'rs');
set(he,'LineWidth',2,'MarkerSize',10,'MarkerFaceColor','r');

% To prevent nearby values from superimposing, we can jitter the data along
% the x-axis. If we don't jitter too much, it will still be clear to which
% group the values belong.
jFactor = (xVals(2) - xVals(1)) / 20;
plot((ones(1,nMale) .* xVals(2)) + (randn(1,nMale).*jFactor),ds.bpf(male),'ko');
plot((ones(1,nFemale) .* xVals(1)) + (randn(1,nFemale).*jFactor),ds.bpf(~male),'ko');

%% Compare bpf of men vs. women using a 2-sample t-test
[h,p,ci,tStats] = ttest2(ds.bpf(male),ds.bpf(~male));

%% Compare bpf of men vs. women using linear regression
const = ones(length(ds),1);
[betaFit,betaCI,resid,residInt,rStats] = regress(ds.bpf,[const,male]);

%% Plot the regression line

ax = axis;
xData = [ax(1),ax(2)];
yData = betaFit(1) + betaFit(2).*(xData-1); % boxplot shifts to 1, 2
line(xData,yData,'Color','r')

%% How do we interpret our beta coefficients?
% beta0 is the model estimate of bpf when male=0, i.e. mean bfp for female
[mean(ds.bpf(~male)), betaFit(1)]

% beta1 is the unit change in bfp for a unit change in maleness, i.e. the
% mean difference in bfp between males and females:
[mean(ds.bpf(male)) - mean(ds.bpf(~male)), betaFit(2)]

%% What is our p-value for the comparison?
[p, rStats(3)]

%% True beauty of regression:
% Of course, there is another explanatory variable luking in the weeds:
% age. A glance at our scatter plot of fig. 1 reveals what appears to be a
% strong effect of age and makes us suspicious that the value for male bpf
% might be driven by the incidental fact that they were a younger sample.

[betaFit2,betaCI2,resid2,residInt2,rStats2] = regress(ds.bpf,[const,ds.age,male]);

% look at the fits and their 95% CIs
[betaFit2, betaCI2]
% an even clearer lack of effect of sex, once we account for age

% plot this regression model:
figure(1);
ax = axis;
xData = [ax(1),ax(2)];
yData = betaFit2(1) + betaFit2(2).*(xData);
line(xData,yData,'Color','r')

% plot model as if m/f difference mattered:
mData = betaFit2(1) + betaFit2(2).*(xData) + betaFit2(3).*[1,1];
fData = betaFit2(1) + betaFit2(2).*(xData) + betaFit2(3).*[0,0];
line(xData,mData,'Color','k','LineStyle','--')
line(xData,fData,'Color','k','LineStyle','-.')

%% Model #3: allow for different slopes for m/f

[betaFit3,betaCI3,resid3,residInt3,rStats3] = regress(ds.bpf,[const,ds.age,male,ds.age.*male]);

[betaFit3, betaCI3]