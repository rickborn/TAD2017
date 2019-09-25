% regVttest.m
%
% Getting the same answer with regression and a t-test
%
% RTB wrote it, 22 Sept. 2017 (hurricane Jose rainy equinox)

% Using data from one of Brian Healy's practicums: bpf.xls
% (originally, 'bpf.dta'; I converted it to an excel file)
%
% Each row is one patient
% Each column is a variable:
% col 1: age in years
% col 2: bpf ( = 'brain parenchymal fraction') How much of your head is brain?
% col 3: gender (M,F)

% Scenario: In a sample of patients with multiple sclerosis (MS), we have
% used a head C/T scan to measure the so-called 'brain parenchymal
% fraction' or 'bpf', which is a commonly used measure of whole brain
% atrophy in MS patients. For each patient, we also have the 'age' and
% 'gender'.
% 
% Scientific Question: Does 'bpf' differ in MS patients who are male vs.
% female?

%% Read in data
ds = readtable('bpf.xls');

% create an indicator variable for gender: 1=male, 0=female
ds.male = strcmp(ds.gender,'M');
nMale = sum(ds.male);
nFemale = sum(~ds.male);

%% Make a boxplot of the bpf broken down by gender

figure
subplot(3,1,1);
boxplot(ds.bpf,ds.male);
hold on

ylabel('brain parenchymal fraction');
set(gca,'XTickLabel',['F';'M'])

% Add mean +/- sem
semMale = std(ds.bpf(ds.male)) / sqrt(nMale);
semFemale = std(ds.bpf(~ds.male)) / sqrt(nFemale);
he = errorbar(xVals,[mean(ds.bpf(~ds.male)),mean(ds.bpf(ds.male))],[semFemale,semMale],'rs');
set(he,'LineWidth',2,'MarkerSize',10,'MarkerFaceColor','r');

% Superimpose raw data on the box plots
% To do this, we need to know the numeric values associated with the two
% groups. This can be obtained using the 'gca' command
xVals = get(gca,'XTick');

% To prevent nearby values from superimposing, we can jitter the data along
% the x-axis. If we don't jitter too much, it will still be clear to which
% group the values belong.
jFactor = (xVals(2) - xVals(1)) / 20;
plot((ones(1,nMale) .* xVals(2)) + (randn(1,nMale).*jFactor),ds.bpf(ds.male),'ko');
plot((ones(1,nFemale) .* xVals(1)) + (randn(1,nFemale).*jFactor),ds.bpf(~ds.male),'ko');

%% Compare bpf of men vs. women using a 2-sample t-test

[h,p,ci,tStats] = ttest2(ds.bpf(ds.male),ds.bpf(~ds.male));

%% Compare bpf of men vs. women using linear regression

const = ones(length(ds.bpf),1);
[betaFit,betaCI,resid,residInt,rStats] = regress(ds.bpf,[const,ds.male]);

modelSpec = 'bpf ~ male';
mdl1 = fitglm(ds,modelSpec,'Distribution','normal');

%% Plot the regression line

ax = axis;
xData = [ax(1),ax(2)];
yData = betaFit(1) + betaFit(2).*(xData-1); % boxplot shifts to 1, 2
line(xData,yData,'Color','r')
title('Model #1: male');

%% How do we interpret our beta coefficients?

% beta0 is the model estimate of bpf when male=0, i.e. mean bfp for female
[mean(ds.bpf(~ds.male)), betaFit(1)];

% beta1 is the unit change in bfp for a unit change in maleness, i.e. the
% mean difference in bfp between males and females:
[mean(ds.bpf(ds.male)) - mean(ds.bpf(~ds.male)), betaFit(2)];

%% What is our p-value for the comparison?

[p, rStats(3)]

%% True beauty of regression: lurking co-variates

% Of course, there is another explanatory variable luking in the weeds:
% age. A glance at our scatter plot of fig. 1 reveals what appears to be a
% strong effect of age and makes us suspicious that the value for male bpf
% might be driven by the incidental fact that they were a younger sample.

subplot(3,1,2);
gscatter(ds.age,ds.bpf,ds.gender,'br','xo');
hold on
xlabel('age (yrs)'); ylabel('brain parenchymal fraction');
title('Model #2: age + male');

[betaFit2,betaCI2,resid2,residInt2,rStats2] = regress(ds.bpf,[const,ds.age,ds.male]);

% look at the fits and their 95% CIs
% [betaFit2, betaCI2]
% an even clearer lack of effect of sex, once we account for age

%% Plot this regression model:

ax = axis;
xData = [ax(1),ax(2)];
yData = betaFit2(1) + betaFit2(2).*(xData);
% line(xData,yData,'Color','k')

% same thing using 'fitglm':
modelSpec = 'bpf ~ age + male';
mdl2 = fitglm(ds,modelSpec,'Distribution','normal');

% plot model as if m/f difference mattered:
mData = betaFit2(1) + betaFit2(2).*(xData) + betaFit2(3).*[1,1];
fData = betaFit2(1) + betaFit2(2).*(xData) + betaFit2(3).*[0,0];
line(xData,mData,'Color','b','LineStyle','-')
line(xData,fData,'Color','r','LineStyle','-')

%% Model #3: allow for different slopes for m/f

subplot(3,1,3);
gscatter(ds.age,ds.bpf,ds.gender,'br','xo');
hold on
xlabel('age (yrs)'); ylabel('brain parenchymal fraction');
title('Model #3: age + male + age:male');

[betaFit3,betaCI3,resid3,residInt3,rStats3] = regress(ds.bpf,[const,ds.age,ds.male,ds.age.*ds.male]);

% plot model as if m/f difference mattered:
mData = betaFit3(1) + betaFit3(2).*(xData) + betaFit3(3).*[1,1] + betaFit3(4).*xData;
fData = betaFit3(1) + betaFit3(2).*(xData) + betaFit3(3).*[0,0];
line(xData,mData,'Color','b','LineStyle','-')
line(xData,fData,'Color','r','LineStyle','-')

% using fitglm:
% The formula for the model is expressed in Wilkinson notation:
modelSpec = 'bpf ~ age + male + age:male';
mdl3 = fitglm(ds,modelSpec,'Distribution','normal');

% [betaFit3, betaCI3]