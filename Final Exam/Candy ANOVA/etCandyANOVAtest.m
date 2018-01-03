%% Candy production ANOVA
% This exercise is based on data from the federal reserve, accessed via
% kaggle: https://www.kaggle.com/rtatman/us-candy-production-by-month
%
% Exercise written by RAS 10/17
% Simplified to use data that has been pre-sorted by season: RTB 24 Oct
% 2017

% Questions: How does candy production vary throughout the year?
%            Which months have the greatest candy production?
%            Which season has the greatest candy production?
%% Import and plot data

% NOTE: RTB filled in the 2 NaN's with column means and converted it to an
% excel spreadsheet so the students would not spend time having to figure
% this out.
fileName='candy-production-by-season.xlsx';
ds = dataset('xlsfile',fileName);

seasonsums = [ds.Winter, ds.Spring, ds.Summer, ds.Fall];
seasonMeans=nanmean(seasonsums);
seasonCI=1.96.*(nanstd(seasonsums)/sqrt(length(seasonsums(:,1))));

figure
x=ones(1,length(seasonsums(:,1)));
x1=x+rand(size(x))-0.5;
x2=x+rand(size(x))+.5;
x3=x+rand(size(x))+1.5;
x4=x+rand(size(x))+2.5;
plot([x1',x2',x3',x4'],seasonsums,'o')
xlabel('Season')
ylabel('Candy production')
hold on
% RTB: removing 'CapSize' again.
errorbar(seasonMeans,seasonCI,'.k','LineWidth', 2,'MarkerSize', 20)
legend('Winter', 'Spring','Summer','Fall','Location','south','Orientation','horizontal')

%Is there a visible pattern of candy production by season? Are there any
%other sensible ways to split up the data other than by season?

%% Does candy production differ by season? ANOVA approach

% One way to answer this question is by using One-way ANOVA
% What is our null hypothesis? - all seasons are the same 

[pVal,tbl]=anova1(seasonsums);
% What does the anova table suggest about our data? Are all seasons drawn
% from the same distribution?
obsF=tbl{2,5};

%% Does candy production differ by season? Permutation test

nPerm = 1000;
% First, let's create our population of season values, ignoring 2017
allseasonsums=[seasonsums(1:end-1,1);seasonsums(1:end-1,2);seasonsums(1:end-1,3);seasonsums(1:end-1,4)];

% Then, sample without replacement from these values, separate into fake
% seasons and perform the ANOVA again to get a distribution of F under H0.
seasonTable=ones(45,4);
bootFs=zeros(1,length(nPerm));
for k=1:nPerm
    bootPop=datasample(allseasonsums,length(allseasonsums), 'Replace', false);
    seasonTable=reshape(bootPop,[45,4]);
    [~,bootTbl]=anova1(seasonTable,[],'off');
    bootFs(k)=bootTbl{2,5};
end

%% Plot a distribution of the Fs we got under H0 

figure
histogram(bootFs);
bsAxis = axis;
hold on

% and plot the F stat we observed from our data. 
% line([obsF,obsF],[bsAxis(3),bsAxis(4)],'Color','k');

%% Calculate a 95% CI and plot it

myAlpha=0.05;
bootFsSorted = sort(bootFs);

idxLo = ceil((myAlpha/2) * nPerm);   % index corresponding to lower bound
idxHi = nPerm - idxLo;               % index corresponding to upper bound

F95CIpercentileLow = bootFsSorted(idxLo);
F95CIpercentileHi = bootFsSorted(idxHi);
line([F95CIpercentileLow,F95CIpercentileLow],[bsAxis(3),bsAxis(4)],'Color','r');
line([F95CIpercentileHi,F95CIpercentileHi],[bsAxis(3),bsAxis(4)],'Color','r');

% What is the 95% confidence interval for F statistics under H0?
% Does candy production depend on season? 

%% Testing assumptions
% 1. Equal variances
[p,stats] = vartestn(seasonsums);
% Our different groups have similar variances 

% 2. Normality of the residuals. What this is really trying to get at is 
% that the distribution of Y|X is normal. In our data, we have a 
% categorical predictor variable and a continuous dependent variable. 
% Here residuals will have the same distribution as Y values in each group.
% We can look at each group individually and assess each for normality.
figure
subplot(2,2,1)
h=qqplot(seasonsums(:,1));
% RTB: title names were in double-quotes; I changed to single
title('Winter')
subplot(2,2,2)
h=qqplot(seasonsums(:,2));
title('Spring')
subplot(2,2,3)
h=qqplot(seasonsums(:,3));
title('Summer')
subplot(2,2,4)
h=qqplot(seasonsums(:,4));
title('Fall')

% Is Y|X normally distributed for all seasons? If not, what can we do about
% it?
%% Alternatively, we can look at all residuals together
winter=seasonsums(:,1);
spring=seasonsums(:,2);
summer=seasonsums(:,3);
fall=seasonsums(:,4);
winterResid=winter-nanmean(winter);
springResid=spring-nanmean(spring);
summerResid=summer-nanmean(summer);
fallResid=fall-nanmean(fall);
allResid=[winterResid;springResid;summerResid;fallResid];
figure
subplot(4,4,1)
qqplot(allResid)
title('Untransformed')
%% Transform our data (all residuals)
% For box cox transformations, our data have to all be positive
allResid=allResid+abs(min(allResid))+0.1.*abs(min(allResid));
allResid(isnan(allResid(:))) = [];
[allResidOpt,lambda] = boxcox(allResid);
subplot(4,4,2)
qqplot(allResidOpt)
title(['Optimum Box Cox lambda=' num2str(lambda)]);
lambdas=[-3:0.5:3];
for y=1:13
    [allResidy] = boxcox(lambdas(y),allResid);
    subplot(4,4,y+2)
    qqplot(allResidy)
    title(['lambda=' num2str(lambdas(y))])
end
    
%% Transform our data (by season):

seasonsumsT=1./(seasonsums);
seasonsumsT=exp(seasonsums);

figure
subplot(2,2,1)
h=qqplot(seasonsumsT(:,1));
title('Winter')
subplot(2,2,2)
h=qqplot(seasonsumsT(:,2));
title('Spring')
subplot(2,2,3)
h=qqplot(seasonsumsT(:,3));
title('Summer')
subplot(2,2,4)
h=qqplot(seasonsumsT(:,4));
title('Fall')