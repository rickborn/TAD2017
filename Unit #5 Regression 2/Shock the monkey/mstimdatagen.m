function [rawData] = mstimdatagen(fn,n)

% mstimdatagen: generate trial-based data from proportion data
% 
% This function reads in an excel file containing proportion data from
% microstim experiments described in Salzman et al. 1992, fig. 4, and
% generates the appropriate number of yes/no trials
%
% Inputs:
% - fn: file name, e.g. es5b_prop.xlsx
% - n: number of trials on which each data point is based
%
% Ouputs:
% - rawData: rows for each trial
%
% RTB wrote it, 12 September 2016

[num] = xlsread(fn);
[nRows,nCols] = size(num);

% init array for raw data
rawData = zeros(nRows * n, nCols);

for thisRow = 1:nRows
    nPDchoices = round(n * num(thisRow,3));
    thisChoice = [ones(nPDchoices,1);zeros(n-nPDchoices,1)];
    thisLabel = repmat(num(thisRow,1:2),n,1);
    thisBlock = [thisLabel, thisChoice];
    startRow = ((thisRow - 1) * n) + 1;
    stopRow = startRow + (n - 1);
    rawData(startRow:stopRow,:) = thisBlock;
end

% now scramble the row order to simulate real experiment
newInd = randperm(length(rawData));
rawData = rawData(newInd,:);

end

% To do logistic regression using the glm:
% [b,dev,stats] = glmfit(rawData(:,[1,2]),rawData(:,3),'binomial','link','logit');