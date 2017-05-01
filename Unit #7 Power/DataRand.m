function [SimulatedExpData1,SimulatedExpData2] = DataRand()

% DataRand: generates data that has the same mean and standard deviation as
% the original prostate cancer data set

% This is the real data set
load prostatecancerexpdata

nReps = size(independentData,2); %This is the number of replicates for every probe
xMu = mean(independentData,2);   %We compute the mean expression of the replicates
xStd = std(independentData')';   %We compute the standard deviation for each gene in the independent data set  

% Now we have two parameters (mean and standard deviation)to specify a new
% distribution from which new gene expression data will be generated

% normrnd takes 2 parameters and can take an array of means Xmu and and array
% of standard deviations Xstd to create a resulting array that has the same
% dimensions. Note that the array of means Xmu and the array of standarad
% deviations Xstd have to have the same number of elements

SimulatedExpData1 = normrnd(repmat(xMu,[1 nReps]), repmat(xStd,[1 nReps]));
SimulatedExpData2 = normrnd(repmat(xMu,[1 nReps]), repmat(xStd,[1 nReps]));

