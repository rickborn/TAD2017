% tSNEdemo.m: a kind of failure mode of t-SNE
%
% RTB wrote it, 25 July 2018, home with Bentley after his biopsy

% MATLAB topics covered:
% 1. using 'mvnrnd' to generate simulated data from bivariate normal distributions
% 2. plotting labeled data using 'gscatter'
% 3. using 'tsne' for dimensionality reduction

% tSNE is short for: t-distributed stochastic neighbor embedding
% See the wikipedia page for a brief explanation of how it works:
% https://en.wikipedia.org/wiki/T-distributed_stochastic_neighbor_embedding

% This is demo is based off of a StackExchange answer by Erich Schubert:
% https://stats.stackexchange.com/questions/263539/clustering-on-the-output-of-t-sne/264647#264647

%% Load some bivariate normal data:
load labeledData

% my default font size:
% myFS = 'default';
myFS = 12;

%% Generate new simulated data using 'mvnrnd'

genFlag = 0;
% NOTE: You don't need to do this step if you have already loaded the
% 'labeledData' data set.

if genFlag
    % population means for bivariate normal data
    muRed = [-2,0];
    muBlue = [2,0];
    
    % population covariance matrices
    sigmaRed = [1,0;0,1];
    sigmaBlue = sigmaRed;
    
    % generate random samples from the correpsonding distributions:
    redData = mvnrnd(muRed,sigmaRed,250);
    blueData = mvnrnd(muBlue,sigmaBlue,750);
    
    % concatenate the data and generate labels:
    data = [redData;blueData];
    labels = [repmat({'red'},250,1);repmat({'blue'},750,1)];
    nLabels = [ones(250,1);ones(750,1).*2];
end

%% Plot the data

figure
subplot(2,2,1);
gscatter(data(:,1),data(:,2),labels,'rb','+x',[],'off');
xlabel('X1');
ylabel('X2');
title('Two bivariate normal distributions');
set(gca,'FontSize',myFS);

%% Now do tSNE using default settings

rng default     % for consistency
Y = tsne(data);

subplot(2,2,2);
gscatter(Y(:,1),Y(:,2),labels,'rb','+x',[],'off');
xlabel('tSNE1');
ylabel('tSNE2');
title('tSNE plot with default settings');
set(gca,'FontSize',myFS);

% Run the tsne algorithm several different times, each time plotting the
% results. (Do NOT run 'rng default' before each tsne run!). What do you
% notice? What does the 'S' in 'tSNE' stand for?

%% Now use a 'Perplexity' value that is too low

% Perplexity is the effective number of local neighbors of each point,
% specified as a positive scalar. The default is 30
% 
% Larger perplexity causes tsne to use more points as nearest neighbors.
% Use a larger value of Perplexity for a large dataset. Typical Perplexity
% values are from 5 to 50. In the Barnes-Hut algorithm, tsne uses
% min(3*Perplexity,N-1) as the number of nearest neighbors.

myPerplexity = 10;
rng default
Y = tsne(data,'Perplexity',myPerplexity);

subplot(2,2,3);
gscatter(Y(:,1),Y(:,2),labels,'rb','+x',[],'off');
xlabel('tSNE1');
ylabel('tSNE2');
title(['tSNE plot with Perplexity = ' num2str(myPerplexity) '. Fake clusters!']);
set(gca,'FontSize',myFS);

%% Run with a more appropriate perplexity

myPerplexity = 60;
rng default
Y = tsne(data,'Perplexity',myPerplexity);

subplot(2,2,4);
gscatter(Y(:,1),Y(:,2),labels,'rb','+x',[],'off');
xlabel('tSNE1');
ylabel('tSNE2');
title(['tSNE plot with Perplexity = ' num2str(myPerplexity) '.']);
set(gca,'FontSize',myFS);

% From the Schubert post:
% "Now this is visually pleasing, but not better for analysis. A human
% annotator could likely select a cut and get a decent result; k-means
% however will fail even in this very very easy scenario! You can already
% see that density information is lost: all data seems to live in area of
% almost the same density. If we would instead further increase the
% perplexity, the uniformity would increase, and the separation would
% reduce again.
% 
% In conclusions, use t-SNE for visualization (and try different parameters
% to get something visually pleasing!), but rather do not run clustering
% afterwards, in particular do not use distance- or density based
% algorithms, as this information was intentionally (!) lost.
% Neighborhood-graph based approaches may be fine, but then you don't need
% to first run t-SNE beforehand, just use the neighbors immediately
% (because t-SNE tries to keep this nn-graph largely intact)."

%% New topic: tSNE with different distance metrics

load fisheriris

rng('default') % for reproducibility
Y = tsne(meas,'Algorithm','exact','Distance','mahalanobis');

figure('Name','tSNE with different distance metrics');
subplot(2,2,1)
gscatter(Y(:,1),Y(:,2),species)
title('Mahalanobis')

rng('default') % for fair comparison
Y = tsne(meas,'Algorithm','exact','Distance','cosine');
subplot(2,2,2)
gscatter(Y(:,1),Y(:,2),species)
title('Cosine')

rng('default') % for fair comparison
Y = tsne(meas,'Algorithm','exact','Distance','chebychev');
subplot(2,2,3)
gscatter(Y(:,1),Y(:,2),species)
title('Chebychev')

rng('default') % for fair comparison
Y = tsne(meas,'Algorithm','exact','Distance','euclidean');
subplot(2,2,4)
gscatter(Y(:,1),Y(:,2),species)
title('Euclidean')

%% Other clustering methods: K-means

X = data;
opts = statset('Display','off');
[idx,C] = kmeans(X,2,'Distance','cityblock',...
    'Replicates',5,'Options',opts);

% How accurate were we?
% We originally used a cell array for word labels. But to compare with the
% output of k-means, we want number labels:
% ground truth labels:
nLabels = [ones(250,1);ones(750,1).*2];
accKmeans = sum(nLabels == idx) / length(idx);
% The label assignment is arbitrary and might not match ours. We want the
% first set to be '1'
if accKmeans < 0.5  % need to swap 1s and 2s
    idx(idx==1) = 3;
    idx(idx==2) = 1;
    idx(idx==3) = 2;
    accKmeans = sum(nLabels == idx) / length(idx);
end

figure;
plot(X(idx==1,1),X(idx==1,2),'ro','MarkerSize',12)
hold on
plot(X(idx==2,1),X(idx==2,2),'bo','MarkerSize',12)
plot(C(:,1),C(:,2),'kx',...
     'MarkerSize',15,'LineWidth',3) 
legend('Cluster 1','Cluster 2','Centroids',...
       'Location','NW')
title 'K-means: Cluster Assignments and Centroids'
xlabel('X1');
ylabel('X2');

% superimpose ground truth data
gscatter(X(:,1),X(:,2),labels,'rb','..',[12,12],'off');

tStr = sprintf('Accuracy = %0.1f%%',accKmeans*100);
text(-4,-3,tStr);

%% Other methods: agglomerative hierarchical cluster tree

% generate the tree:
% 2nd argument is method for determining distance between clusters:
%   - 'single' uses the distance between the 2 nearest members (default)
%   - 'complete' uses the distance between the 2 farthest members
%   - 'average' uses the average of all pairs in any 2 clusters
%   - 'centroid' uses the distance between the centroids
%
% According to Hastie and Tibshirani, 'complete' and 'average' perform
% best.

% 3rd argument is the metric for measuring distance. The default is
% 'euclidean', but there are many other possible metrics.
Z = linkage(X,'average','euclidean');

% view the tree
figure
dendrogram(Z);

% cluster the data
% 'c' contains the label for each data point
c = cluster(Z,'maxclust',2);
dendrogram(Z,'ColorThreshold',4);
hold on
ax = axis;
h1 = line([ax(1),ax(2)],[4,4]);
set(h1,'LineStyle','--');
title('Agglomerative hierarchical cluster tree');

% plot the labeled data
figure
gscatter(X(:,1),X(:,2),c,'rb','oo',[12,12],'off');
hold on
% superimpose ground truth labels
gscatter(X(:,1),X(:,2),labels,'rb','..',[12,12],'off');
xlabel('X1');
ylabel('X2');
title('Agglomerative hierarchical cluster tree');

% How accurate were we?
accTree = sum(nLabels == c) / length(c);
tStr = sprintf('Accuracy = %0.1f%%',accTree*100);
text(-4,-3,tStr);