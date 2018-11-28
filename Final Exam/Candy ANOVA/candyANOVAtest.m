%% candyANOVAtest.m: ANOVA bootstrap for final exam
%
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

fileName='candy-production-by-season.xlsx';
ds = readtable(fileName);

% format data for ANOVA
seasonSums = [ds.Winter, ds.Spring, ds.Summer, ds.Fall];
seasonMeans = mean(seasonSums);
seasonCI=1.96.*(std(seasonSums)/sqrt(length(seasonSums(:,1))));

figure
% Jitter the data so that we can see it all
x=ones(1,length(seasonSums(:,1)));
x1=x+rand(size(x))-0.5;
x2=x+rand(size(x))+.5;
x3=x+rand(size(x))+1.5;
x4=x+rand(size(x))+2.5;
% Plot it
plot([x1',x2',x3',x4'],seasonSums,'o')
xlabel('Season')
ylabel('Candy production')
hold on

% RTB: removing 'CapSize' again.
errorbar(seasonMeans,seasonCI,'.k','LineWidth', 2,'MarkerSize', 20)
legend('Winter', 'Spring','Summer','Fall','Location','south','Orientation','horizontal')

%Is there a visible pattern of candy production by season? Are there any
%other sensible ways to split up the data other than by season?

%% Does candy production differ by season? ANOVA approach

% One way to answer this question is by using One-way ANOVA on the data in
% 'seasonSums'

!!! Your code here:
[pVal,tbl]=anova1();

% What does the anova table suggest about our data? Are all seasons drawn
% from the same distribution?
obsF=tbl{2,5};

%% Does candy production differ by season? Permutation test

nPerm = 1000;

% HINT: You want to shuffle all of the data as if you had misplaced all of
% the season labels, but then reformat the data so that it looks like your
% original data set (i.e. separated into four columns of 'fake' seasons),
% then run 'anova1' to get the F-statistic

% variable to hold each permuted F-statistic:
permFs = zeros(nPerm,1);

!!! Your code here:

for k=1:nPerm
    % shuffle the data:
    
    % re-format into original shape:
    
    % perform a 1-way ANOVA:
    % NOTE: You can suppress the display of the ANOVA table:
    % [] = anova1(data,[],'off');
    
    % extract the F-statistic:
    permFs(k) = Tbl{2,5};
end

%% Plot a distribution of the Fs we got under H0 

!!! Your code here

%% Compare your observed F-statistic to the permuted distribution

% What is your p-value?

!!! Your code here?

