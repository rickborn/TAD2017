%% Candy production ANOVA
% This exercise is based on data from the federal reserve, accessed via
% kaggle: https://www.kaggle.com/rtatman/us-candy-production-by-month
%
% Exercise written by RAS 10/17
%
% Questions: How does candy production vary throughout the year?
%            Which months have the greatest candy production?
%            Which season has the greatest candy production?
%% Import and sort data
fileName='candy_production.csv';
ds = dataset('File',fileName,'Delimiter',',');

% It would be nice to have numbers to represent months and seasons 
% rather than dates
for k=1:length(ds.observation_date)
    currentEntry=ds.observation_date{k};
    ds.year(k)=str2double(currentEntry(1:4));
    ds.month(k)=str2double(currentEntry(6:7));
    if ds.month(k)== 12|ds.month(k)== 1|ds.month(k)== 2
        ds.season(k)={'Winter'};
    elseif ds.month(k)>=3 && ds.month(k)<=5
        ds.season(k)={'Spring'};
    elseif ds.month(k)>=6 && ds.month(k)<=8
        ds.season(k)={'Summer'};
    elseif ds.month(k)>=9 && ds.month(k)<=11
        ds.season(k)={'Fall'};
    end
end
    
%% Plot candy production for each month
monthmeans=zeros(1,12);
monthCI=zeros(1,12);
for k=1:12
    monthmeans(k)=mean(ds.IPG3113N(ds.month==k));
    monthCI(k)=1.96.*std(ds.IPG3113N(ds.month==k))./sqrt(length(ds.IPG3113N(ds.month==k)));
end

figure
errorbar(monthmeans,monthCI,'.','MarkerSize',30,'Color','k','CapSize',20,'LineWidth',2)
xlabel('Month')
ylabel('Candy Production')

% Which months have the highest candy production? Do you have any idea why
% that might be?
%% Plot candy production by year
% Let's get a list of candy production broken up by year 
yearlist=unique(ds.year);
yearsums=zeros(1,length(yearlist));
for k=1:length(yearlist)
    yearsums(k)=sum(ds.IPG3113N(double(ds.year)==yearlist(k)));
end

%Because we don't have all of 2017, let's not plot it
yearsums=yearsums(1:end-1);

yearMeans=mean(yearsums);
yearCI=1.96.*std(yearsums)/sqrt(length(yearsums(:,1)));

figure
plot(yearlist(1:end-1)',yearsums,'ok')
xlabel('Year')
ylabel('Candy production')

% Describe the data. How has US candy production varied over the past 30+
% years? Do you see any trends?
%% Plot candy production by season approach #1: by calendar year
% This approach puts December of a given year with the previous January and
% February
yearlist=unique(ds.year);
seasons={'Winter', 'Spring','Summer','Fall'};
seasonsums=zeros(length(yearlist),length(seasons));
%This code will code seasons by year. This has the disadvantage of grouping
%January and February with the following December. Maybe it would make more
%sense to group last year's December with January and February
for k=1:length(yearlist)
    for j=1:length(seasons)
        seasonsums(k,j)=sum(ds.IPG3113N(double(ds.year==yearlist(k)) + double(ismember(ds.season,seasons{j}))==2));
    end
end
%You may notice we have two 0s in our dataset at this point. This is
%because we don't have the full totals for this fall and this winter. Let's
%remove those 0s now. 
seasonsums(46,[1,4])=NaN;

%% Alternate approach(Dec-Jan-Feb together):
%This approach keeps a single winter together by grouping December from the
% previous year with January and February of the next year

%This approach relies on a dataset that is sorted in a particular way
dsSorted=sortrows(ds,[3,4],'ascend');
seasonsums=NaN(length(dsSorted.year),length(seasons));
for k=1:length(dsSorted.year)-2
    for j=1:length(seasons)
        if strcmp(ds.season{k},seasons{j}) && strcmp(ds.season{k+1},seasons{j}) && strcmp(ds.season{k+2},seasons{j})
            seasonsums(k,j)=sum(ds.IPG3113N(k:k+2));
        end
    end
end

%Remove the NaNs in each columns and then add back one NaN to the end of
%winter and fall to account for us missing fall and winter 2017 data
winter=seasonsums(:,1);
winter(isnan(winter(:,1)),:) = [];
winter=[winter;NaN];
spring=seasonsums(:,2);
spring(isnan(spring(:,1)),:) = [];
summer=seasonsums(:,3);
summer(isnan(summer(:,1)),:) = [];
fall=seasonsums(:,4);
fall(isnan(fall(:,1)),:) = [];
fall=[fall;NaN];

%combine these columns into the new season matrix
seasonsums=[winter,spring,summer,fall]
%% Plot season means and CI
% Upload a plot of candy production by season, including both the raw data
% and the CI
seasonMeans=nanmean(seasonsums);
seasonCI=1.96.*nanstd(seasonsums)/sqrt(length(seasonsums(:,1)));

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
errorbar(seasonMeans,seasonCI,'.k','LineWidth', 2,'CapSize', 20,'MarkerSize', 20)
legend('Winter', 'Spring','Summer','Fall','Location','south','Orientation','horizontal')

%Is there a visible pattern of candy production by season? Are there any
%other sensible ways to split up the data other than by season?
%% Does candy production differ by season?
% One way to answer this question is by using One-way ANOVA
% What is our null hypothesis? - all seasons are the same 

[p,tbl]=anova1(seasonsums);
% What does the anova table suggest about our data? Are all seasons drawn
% from the same distribution?
obsF=tbl{2,5};


% Let's do this with bootstrapping
nBoot=1000;
% First, let's create our population of season values, ignoring 2017
allseasonsums=[seasonsums(1:end-1,1);seasonsums(1:end-1,2);seasonsums(1:end-1,3);seasonsums(1:end-1,4)];

% Then, sample with replacement from these values, separate into fake
% seasons and perform the ANOVA again to get a distribution.
bootPop=ones(length(allseasonsums));
seasonTable=ones(45,4);
bootFs=zeros(1,length(nBoot));
for k=1:nBoot
    bootPop=datasample(allseasonsums,length(allseasonsums), 'Replace', false);
    seasonTable=reshape(bootPop,[45,4]);
    [bootp,bootTbl]=anova1(seasonTable,[],'off');
    bootFs(k)=bootTbl{2,5};
end

%% Plot a distribution of the Fs we would get under H0 
% and plot the F stat we observed from our data. 
figure
histogram(bootFs);
bsAxis = axis;
hold on
line([obsF,obsF],[bsAxis(3),bsAxis(4)],'Color','g');

%Calculate a 95% CI and plot that
myAlpha=0.05
bootFsSorted = sort(bootFs);

idxLo = ceil((myAlpha/2) * nBoot);   % index corresponding to lower bound
idxHi = nBoot - idxLo;               % index corresponding to upper bound

F95CIpercentileLow = bootFsSorted(idxLo);
F95CIpercentileHi = bootFsSorted(idxHi);
line([F95CIpercentileLow,F95CIpercentileLow],[bsAxis(3),bsAxis(4)],'Color','r');
line([F95CIpercentileHi,F95CIpercentileHi],[bsAxis(3),bsAxis(4)],'Color','r');

% What is the 95% confidence interval for F statistics under H0?
% Does candy production depend on season? 
%% Post-hoc testing?

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
h=qqplot(seasonsums(:,1))
title("Winter")
subplot(2,2,2)
h=qqplot(seasonsums(:,2))
title("Spring")
subplot(2,2,3)
h=qqplot(seasonsums(:,3))
title("Summer")
subplot(2,2,4)
h=qqplot(seasonsums(:,4))
title("Fall")

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
title("Optimum Box Cox lambda="+num2str(lambda));
lambdas=[-3:0.5:3];
for y=1:13
    [allResidy] = boxcox(lambdas(y),allResid);
    subplot(4,4,y+2)
    qqplot(allResidy)
    title("lambda="+num2str(lambdas(y)))
end
    
%% Transform our data (by season):
seasonsumsT=1./(seasonsums)
seasonsumsT=exp(seasonsums)

figure
subplot(2,2,1)
h=qqplot(seasonsumsT(:,1))
title("Winter")
subplot(2,2,2)
h=qqplot(seasonsumsT(:,2))
title("Spring")
subplot(2,2,3)
h=qqplot(seasonsumsT(:,3))
title("Summer")
subplot(2,2,4)
h=qqplot(seasonsumsT(:,4))
title("Fall")