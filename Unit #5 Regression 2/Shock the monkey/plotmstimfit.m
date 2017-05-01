function [equivVisualStimulus] = plotmstimfit(fn,b)

% plotmstimfit: Uses parameters from logistic regression to plot fits
%
% EV = plotmstimfit(fn,b)
%
% Inputs:
% - fn: file name containing proportions of PD choices: e.g. js21aProp.xlsx
% - b: beta parameters for fit
%
% Outputs:
% - EV: equivalent visual stimulus
%
% RTB wrote it 14 September 2016

[num] = xlsread(fn);
[nRows,~] = size(num);

% plot raw data
plot(num(1:nRows/2,2), num(1:nRows/2,3), 'ko');
hold on;
plot(num(nRows/2 + 1:end,2), num(nRows/2 + 1:end,3), 'ro','MarkerFaceColor','r');

% calculate regression lines from fit coefficients
coh = [floor(min(num(:,2))):0.1:ceil(max(num(:,2)))];
const = ones(size(coh));
mstim = zeros(size(coh));
Pnostim = 1 ./ (1 + exp(-(b(1).*const + b(2).*mstim + b(3).*coh)));
plot(coh,Pnostim,'k-')

% allow for possibility of a significant interaction term
if length(b) == 4
    Pstim = 1 ./ (1 + exp(-(b(1).*const + b(2).*const + b(3).*coh + b(4).*coh)));
else
    Pstim = 1 ./ (1 + exp(-(b(1).*const + b(2).*const + b(3).*coh)));
end
plot(coh,Pstim,'r-');
xlabel('Stimulus Strength (%Coh)');
ylabel('Proportion PD');

equivVisualStimulus = b(2) / b(3);  % Why does this work?
eqv = round(equivVisualStimulus * 10) / 10;
tStr = [fn,': Equiv. Visual Stimulus = ',num2str(eqv), '%coh'];
title(tStr);

% NOTE: You could also find the eqv using a brute force approach:
EV = mean(coh(Pnostim < 0.505 & Pnostim > 0.495)) - ...
     mean(coh(Pstim < 0.505 & Pstim > 0.495));

end

