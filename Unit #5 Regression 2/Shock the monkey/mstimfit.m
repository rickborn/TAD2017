function [b,dev,stats] = mstimfit(rawData)

% mstimfit: performs logistic regression on mstim data
% 
% RTB wrote it

% columns: Mstim, Coh, Choice


% First fit full model to see if interaction term is significant
inter = rawData(:,1) .* rawData(:,2);
[b,dev,stats] = glmfit([rawData(:,[1,2]),inter],rawData(:,3),'binomial','link','logit');

if stats.p(4) > 0.05
    [b,dev,stats] = glmfit(rawData(:,[1,2]),rawData(:,3),'binomial','link','logit');
end

end

