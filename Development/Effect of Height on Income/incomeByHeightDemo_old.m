% incomeByHeightDemo.m
%
% RTB trying to emulate Gelman & Nolan demo (pp. 49-51)

%% Load and plot data

% cd 'C:\usr\rick\doc\Committees\PIN\PIN Director\Courses\Stats\TAD Fall 2017\TAD2017\Last class'
fileName = 'IncomeByHgtData.xlsx';
ds = dataset('xlsfile',fileName);

figure

% We want to use actual height for regression, but jittered height for
% plotting:
plot(ds.Hgt, ds.Income, 'b.');
hold on
% draw the least-squares regression line:
lsline
xlabel('Height (inches)');
ylabel('Income (thousands of $)');

%% Do simple regression
modelspec1 = 'Income ~ Hgt';
mdl1 = fitglm(ds,modelspec1,'Distribution','normal');
[b1,dev1,stats1] = glmfit(ds.Hgt,ds.Income);

%% Lurking variable: sex
% men taller than women and make more $$$

ds.Male = logical(ds.Male);
h1 = plot(ds.Hgt(ds.Male),ds.Income(ds.Male),'k+');
h2 = plot(ds.Hgt(~ds.Male),ds.Income(~ds.Male),'ro');
legend([h1,h2],{'Male','Female'},'Location','NorthWest');

%% Add sex to regression
modelspec2 = 'Income ~ Hgt + Male';
mdl2 = fitglm(ds,modelspec2,'Distribution','normal');
[b2,dev2,stats2] = glmfit([ds.Hgt,ds.Male],ds.Income);

%% Show separate regression lines for men and women

% re-plot original data:
figure
h1 = plot(ds.Hgt(ds.Male),ds.Income(ds.Male),'k+');
hold on
h2 = plot(ds.Hgt(~ds.Male),ds.Income(~ds.Male),'ro');
%legend([h1,h2],{'Male','Female'},'Location','NorthWest');
xlabel('Height (inches)');
ylabel('Income (thousands of $)');

%% Get predictions and CIs for men:

xVals = (min(ds.Hgt):max(ds.Hgt))';
men = ones(size(xVals));
const = zeros(size(xVals));

[yMen,yMenCI] = predict(mdl2,[xVals,men]);
% Plot 'em
plot(xVals,yMen,'k-','LineWidth',2);
plot(xVals,yMenCI(:,1),'k--');
plot(xVals,yMenCI(:,2),'k--');

% Get predictions and CIs for women:
[yWomen,yWomenCI] = predict(mdl2,[xVals,~men]);
% Plot 'em
plot(xVals,yWomen,'r-','LineWidth',2);
plot(xVals,yWomenCI(:,1),'r--');
plot(xVals,yWomenCI(:,2),'r--');

legend([h1,h2],{'Male','Female'},'Location','NorthWest');