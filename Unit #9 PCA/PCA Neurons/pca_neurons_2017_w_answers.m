%% PCA on simulated neurons
%
% Guidance is provided via comments and example code below.
%
% To make this script easier to use, the task is broken into  sections,
% each section with a bold header (those starting with '%%') can be run by
% themselves (without running the entire script) using "ctrl + enter"
% (Windows) or "command + enter" (MAC). Just place your cursor within one of
% these sections (the section will become highlighted) to allow this functionality.
% 
% AVB & SLH 4/2016 editted by LND 2017

%% Close figures and clear workspace
clear;      % Delete all  variables in workspace, you will lose unsaved variables
close all;  % Close all of the open figure windows
%% Load data
% This line loads three variables: data, stim, time
load('pca_data.mat')

% data is a 58x5000 matrix, Neurons x Time Points
% Each row of data is the PSTH of a neuron's response to the stimuli
% stim is a 1x5000 column vector of the Stimuli over time
% time is a 1x5000 column vector of the Time in second
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% PART ONE %%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Plot data for five neurons
% This section of the code plots the stimulus and the responses for the
% first five neurons in the same figure.  You can copy and change this code
% to create your other figures.

% If you can't see the data in the figure, maximize the figure so you can. 

figure
ax(1) = subplot(6,1,1); % subplot allows you to plot multiple graphs in the same figure
plot(time,stim','r') % Plot the stimulus in red ('r')
ylabel('Odor concentration')
title('Stimulus')
for i = 1:5 % Loop over first 5 neurons to plot their responses
    ax(i+1) = subplot(6,1,i+1);
    plot(time,data(i,:)) % This plots row i of data
    title(['Response of neuron ',num2str(i)])
end
xlabel('Time (seconds)')
ylabel('Response (spikes/sec)')


% Now plot more example neurons

% Then try to plot the responses of all neurons together in a single plot
figure, waterfall(data);
xlabel('Time (seconds)')
ylabel('Neuron #')
zlabel('Response (spikes/sec)')

figure, imagesc(data);
xlabel('Time (seconds)')
ylabel('Neuron #')
cb = colorbar;
cb.Label.String = 'Response (spikes/sec)';

% Play with different colormaps to see what brings out the structure
colormap('gray')
colormap('hot')
colormap('default')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% PART TWO %%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Explore ways to compare how similar the responses are of neuronal pairs.
% Think about what would be good similarity metrics to compare neuronal
% responses. Try to implement these if you want.


%% Generate the covariance matrix
% Use the 'cov' function to calculate the covariance matrix.  The 'cov'
% function automatically centers the data so you do not need to do this. 
% This function will return a matrix that has variance along the diagonal
% entries and the covariance in the off-diagonal entries.

% Note that for this data we could either calculate a covariance matrix of
% how different neurons covary together or of how different time points
% covary together. (i.e. we calculate either how the rows or columns of the
% data matrix covary). Recall that our goal in the end is to use PCA to
% reduce the dimensionality of the data from 58 neurons to a smaller number
% of ‘principal component’ neurons so think about which covariance matrix
% is most useful to us if we want to achieve this goal.

% Hint: Type 'doc cov' in the command window to look up how the 'cov'
% function works. 

% You have done this correctly if dataCov(1,1) = 1.0474

% Replace the [] with your own code below.

% =======================
% Insert/Modify code here
% By default, looks at covariance between columns, which, in our case is
% time. To get it to do rows, use transpose:
dataCov = cov(data');

% =======================

%% Plot the covariance matrix
% You can use the 'imagesc' function to visualise the covariance matrix.
% Calling the 'colorbar' function adds the color scale on your graph.
% 

% =======================
% Insert/Modify code here
figure, imagesc(dataCov);
cb = colorbar;
cb.Label.String = 'Covariance between neurons';
% =======================


%% Cluster the covariance matrix
% use the kmeans function with the covariance matrix as input. Think about
% how many clusters (k) to input.

% The output will be a list of cluster IDs. Find all the neurons that
% belong to cluster 1, cluster 2, and so on. Plot the responses of all the
% neurons in cluster 1 and compare them to the stimulus vs. time plot. Use
% similar plotting approaches to those used in Part 1.

for nClusters = 3
    % Insert code here
    %nClusters = 2;
    idx = kmeans(dataCov,nClusters);
    
    % sort by k-means index
    dataClustered = [data(idx==1,:); data(idx==2,:); data(idx==3,:); ...
        data(idx==4,:); data(idx==5,:)];
    
    % figure,
    % ax(1) = subplot(2,1,1); % subplot allows you to plot multiple graphs in the same figure
    % plot(time,stim','r') % Plot the stimulus in red ('r')
    % ylabel('Odor concentration')
    % title('Stimulus')
    % subplot(2,1,2)
    % imagesc(dataClustered);
    % xlabel('Time (seconds)')
    % ylabel('Neuron #')
    % cb = colorbar;
    % cb.Label.String = 'Response (spikes/sec)';
    
    % try appending the stimulus row to the sorted data:
    figure
    maxVal = round(max(data(:)));
    D = [stim .* maxVal; dataClustered];
    imagesc(D);
    xlabel('Time (seconds)')
    ylabel('Neuron #')
    cb = colorbar;
    cb.Label.String = 'Response (spikes/sec)';
    title(['Using ',num2str(nClusters), ' clusters'])
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% PART THREE %%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Perform PCA
% Use the 'pca' function to run PCA on the data matrix. Your goal is to
% analyze the relationships among neurons, not the relationships among
% time points. Your hypothesis is that the responses of all 58 neurons
% can be reduced to linear combinations of a few orthogonal basis
% functions, where each basis function can be conceptualized as a
% distinct "response type". You should set up the PCA so that you obtain
% 58 PCs; you hypothesize that only a few of these PCs are needed to
% explain most of the variance in the data.

% The 'pca' function returns several values. For us the important ones are:
% 
%   'score' - This is a matrix containing the representation of the data
%       set in PC space. Essentially, this is the data set after it has
%       been rotated. Recall that the original data set consisted of 58
%       vectors, with each vector representing a neural response measured
%       at 5000 time points; this new matrix therefore also consists of 58
%       vectors measured at 5000 time points. Note that there is some
%       confusion of terminology in the literature.  Sometimes the PC
%       scores are referred to as the "Principal components (PCs)" while
%       other times our new axes (the eigenvectors) are called the
%       "Principal components". We prefer the second usage and hence we
%       will refer to the scores as "PC scores" and to the axes as PCs but
%       you should be familiar with both.
%   'explained' - This is a list of numbers quantifying the percentage of
%       the variance in the data explained by each of the PCs, in 
%       descending order of variance explained.
%   'coeff' - This is a matrix that quantifies the importance of each
%       variable (here, each neuron in the original dataset) in accounting 
%       for the variability of the associated PC. (This matrix gets the 
%       name 'coeff' because it contains the correlation coefficients 
%       between the data in PC space and the original data.) These values 
%       are also known as loadings. Each column of coeff contains 
%       coefficients for one principal component, and the columns are in 
%       descending order of component variance.

% Read the help documentation on pca for further information. 

% If you're confused, we recommend first performing PCA and the plotting
% steps below and then coming back to revise your code based on your
% results. Plotting the outputs from the pca function can be helpful for
% understanding what pca is doing.

% You have done this correctly if coeff(1,1) = 0.0724

% =======================
% Insert/Modify code here

% Same as for 'cov', we're looking for relationships across rows, whereas
% default is across columns. Hence the transpose:
[coeff,score,~,~,explained,~] = pca(data');

% =======================


%% Plot explained variance (~Scree plot)
% Use the output from the 'pca' function above to make a plot of the
% different PC contributions to explained variance in the data.
%
%   hint: To make it easy to see the variance explained by each pc when you
%       plot 'explained' also pass '-o' to the plot function, like this
%       example: plot(explained,'-o')

% =======================
% Insert/Modify code here
figure, plot(explained,'-o');
xlabel('PC #');
ylabel('Variance explained');
% =======================


%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% PART FOUR %%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Use the output of 'pca' to plot the first six principal component scores
% The first six PC scores are the first six vectors in the score matrix,
% where each vector is a list of 5000 numbers.

% =======================
% Insert/Modify code here
score = score'; % flip back to 58 x 5000
figure
ax(1) = subplot(6,1,1); % subplot allows you to plot multiple graphs in the same figure
plot(time,stim','r') % Plot the stimulus in red ('r')
ylabel('Odor concentration')
title('Stimulus')
for i = 1:5 % Loop over first 5 neurons to plot their responses
    ax(i+1) = subplot(6,1,i+1);
    plot(time,score(i,:)) % This plots row i of data
    title(['Response of PC',num2str(i)])
end
xlabel('Time (seconds)')
ylabel('Score')
% =======================


%% Find the covariance matrix of data in the PC space and plot it
% Use the 'imagesc' function and the 'colorbar' function for plotting. You
% should have 58 PCs, so this should be a 58-by-58 matrix (i.e. the same
% size as the previous covariance matrix you plotted).
% 
% You have done this correctly if the first entry in the matrix = 47.1669

% =======================
% Insert/Modify code here
scoreCov = cov(score');
figure, imagesc(scoreCov);
cb = colorbar;
cb.Label.String = 'Covariance between PC Scores';
% =======================


%% Make a 3D plot of each neuron's loadings for the first 3 PCs
% Use 'plot3' to make a 3D plot. Plot each loading as a discrete dot or
% circle for clarity, and please label the axes. Remember the loadings
% correspond to the 'coeff' output of the 'pca' function.
%
%   hint: type 'doc plot3' for info about how to label the axes 
%   hint: pass 'o' to the plot3 function, to plot each loading as a 
%       discrete dot, like this example: 
%       plot3(x_variable, y_variable, z_variable,'o') 
%   hint: using 'grid on' might make your graph
%       more easily viewable

% =======================
% Insert/Modify code here
% coeff is 58 x 58
% Each column of coeff contains coefficients for one principal component, 
% and the columns are in descending order of component variance.
figure, 
plot3(coeff(idx==1,1),coeff(idx==1,2),coeff(idx==1,3),'ro');
hold on
plot3(coeff(idx==2,1),coeff(idx==2,2),coeff(idx==2,3),'go');
plot3(coeff(idx==3,1),coeff(idx==3,2),coeff(idx==3,3),'bo');
xlabel('PC1 loading')
ylabel('PC2 loading')
zlabel('PC3 loading')
grid on


% =======================


%% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% PART FIVE %%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Perform PCA to analyze the relationship between time points
% Now instead of analyzing the relationship between neurons, we are 
% interested in how the population as a whole evolves over time. Your 
% goal is to visualize that population trajectory.

% Make sure you still have the output of your PCA in the workspace or redo
% the PCA


% This code uses a for loop to show the temporal evolution of the stimulus
% in a specified step size. All previous time points are shown in black and
% the most recent past time points are shown in red. One subplot is used
% for the stimulus. Fill in the code to make in the second subplot a 3D
% plot of the first 3 PC scores. At the end of each iteration of the for
% loop, a pause is included to allow you inspect the figure. A key press is
% required to start the next iteration. Play with the step size if you want
% faster or slower steps through the temporal evolution of stimulus and
% neural activity.

stepSize = 50;  % stepping through time
grayVal = 0.7;  % ghost color

figure; hold on;
for i = stepSize+1:stepSize:length(stim)
    ax1 = subplot(2,1,1);
    cla;
    plot(time(1:i),stim(1:i),'-','Color',[grayVal,grayVal,grayVal],'LineWidth',1);
    hold on;
    plot(time(i-stepSize:i),stim(i-stepSize:i),'r-','LineWidth',3);
    axis(ax1,[0 time(end) 0 1])
    xlabel('Time (seconds)');
    ylabel('Stimulus concentration');
    
    ax2 = subplot(2,1,2);
    cla;
    % Insert code to make a 3D plot of the first 3 PCs, like is done in the
    % same loop for stimulus
    plot3(score(1:i,1),score(1:i,2),score(1:i,3),'o','Color',[grayVal,grayVal,grayVal]);
    drawnow
    hold on
    hp=plot3(score(i-stepSize:i,1),score(i-stepSize:i,2),score(i-stepSize:i,3),'ro');
    set(hp,'MarkerFaceColor','r');
    axis(ax2,[-10 20 -10 10 -10 10])
    xlabel('PC1');
    ylabel('PC2');
    zlabel('PC3');
    pause;
end



%% Extension problems
% If you feel confident or would like to gain additional practice, please
% continue by answering the extension problems as outlined in the
% instructions. Some problems will require more coding, you may complete
% this below.
