% incomeByHeightDemo.m
%
% RTB trying to emulate Gelman & Nolan demo

clearall
nSamp = 100;

% From web:
muHgtInF = 64;      % mean height of adult women in inches
varHgtInF = 9;      % variance (sigma = 3)
muHgtInM = 70;      % mean height of adult women in inches
varHgtInM = 12.5;  % sigma

% RTB making up
muIncomeF = 80;     % in thousands
varIncomeF = 60;
muIncomeM = 85;
varIncomeM = 75;

% Need to come up with some plausible co-var numbers for incomve v. height
% NOTE: I can play with the covariances to illustrate various principles.
% For example, values of ~1.2 give an interesting case where the effect of
% height is purely driven by the lurking variable of sex. If I make the
% individual covariances stronger, we can still see an effect when we
% include sex in our regression model.
covHgtIncM = 1.2;
covHgtIncF = 1.2;

muMale = [muHgtInM,muIncomeM];
covMale = [varHgtInM,covHgtIncM;covHgtIncM,varIncomeM];
% generate some fake data for men:
Rmale = mvnrnd(muMale,covMale,nSamp);

% Same for women:
muFemale = [muHgtInF,muIncomeF];
covFemale = [varHgtInF,covHgtIncF;covHgtIncF,varIncomeF];
Rfemale = mvnrnd(muFemale,covFemale,nSamp);

figure
plot(Rmale(:,1),Rmale(:,2),'k+');
hold on
plot(Rfemale(:,1),Rfemale(:,2),'ro');
xlabel('Height (inches)');
ylabel('Income (thousands of $)');

% concatenate the two data sets
R = [Rmale;Rfemale];
G = [ones(nSamp,1);zeros(nSamp,1)];

% do simple regression: 
[b1,dev1,stats1] = glmfit(R(:,1),R(:,2));

% lurking variable may be sex: men taller than women and make more
% add sex to regression
[b2,dev2,stats2] = glmfit([R(:,1),G],R(:,2));

% create a data set by randomizing rows
rowIDs = randperm(nSamp*2)';
D = [R(rowIDs,:),G(rowIDs,:)];

% write data to an excel spreadsheet
% first, let's do some rounding so values look more plausible
% Height to nearest inch:
% D(:,1) = round(D(:,1));
% % Income to nearest $
% D(:,2) = round(D(:,2) .* 1000) ./ 1000;
% xlswrite('IncomeByHgtData.xlsx',D);