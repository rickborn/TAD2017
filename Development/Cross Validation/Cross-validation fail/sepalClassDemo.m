% sepalClassDemo.m
%
% copied from the documentation of MATLAB's 'classify' function

%% Load some data and create a subset

% Fisher's iris data consists of measurements on the sepal length, sepal
% width, petal length, and petal width for 150 iris specimens. There are 50
% specimens from each of three species.
load fisheriris

% To see the names of the 3 species:
unique(species)

% For simplicity, we will only classify 'versicolor' and 'virginica'
sepalLength = meas(51:end,1);
sepalWidth = meas(51:end,2);
group = species(51:end);

%% Plot the data

figure
h1 = gscatter(sepalLength,sepalWidth,group,'rb','v^',[],'off');
set(h1,'LineWidth',2)
legend('Fisher versicolor','Fisher virginica',...
       'Location','NW')
   
%% Classify a grid of measurements on the same scale:

ax = axis;
[X,Y] = meshgrid(linspace(ax(1),ax(2)),linspace(ax(3),ax(4)));
X = X(:); Y = Y(:);
[C,err,P,logp,coeff] = classify([X Y],[sepalLength sepalWidth],...
                                group,'Quadratic');
                            
%% Visualize the classification:

hold on;
gscatter(X,Y,C,'rb','.',1,'off');
K = coeff(1,2).const;
L = coeff(1,2).linear;
Q = coeff(1,2).quadratic;
% Function to compute K + L*v + v'*Q*v for multiple vectors
% v=[x;y]. Accepts x and y as scalars or column vectors.
f = @(x,y) K + [x y]*L + sum(([x y]*Q) .* [x y], 2);

h2 = ezplot(f,ax);
set(h2,'Color','m','LineWidth',2)
axis(ax)
xlabel('Sepal Length (cm)')
ylabel('Sepal Width (cm)')
title('{\bf Classification with Fisher Training Data}')

%% Perform K-fold cross-validation

K = 10;
indices = crossvalind('Kfold',group,K);
% tell the classifier the right answers
cp = classperf(group);
predictors = [sepalLength,sepalWidth];

for k = 1:K
    test = (indices == k);
    train = ~test;
    
    class = classify(predictors(test,:),predictors(train,:),group(train,:),...
                     'quadratic');
    classperf(cp,class,test)
end
cp.ErrorRate

%% Discriminant analysis on the full data set

% First, let's make a scatter plot matrix of all of the predictor
% variables:
figure, plotmatrix(meas);

% This is nice, but unfortunately it does not allow us color the dots
% according to which group they below to. We can do this for some example
% plots, such as the one in row 2, column 1:
figure, gscatter(meas(:,1),meas(:,2),species);

% or the one in row 3 of column 1:
figure, gscatter(meas(:,1),meas(:,3),species);

%% Let's run cross-validation at see how well we do:

K = 10;
indices = crossvalind('Kfold',species,K);
cp = classperf(species);
for k = 1:K
    test = (indices == k); train = ~test;
    class = classify(meas(test,:),meas(train,:),species(train,:),'linear');
    classperf(cp,class,test);
end
cp.ErrorRate