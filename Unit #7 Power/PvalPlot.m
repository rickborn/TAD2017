function [T0,T0sd] = PvalPlot(P,pFlag)

% PVALPLOT.M: Plot of P-values according to Schweder & Spjotvoll 1982
%
% [T0] = PvalPlot(P,pFlag)
%
% Inputs:
% - P: a vector of p-values from a large number of comparisons
% - pFlag: plots a line fit to the lower part of the p-value
%
% Outputs:
% - T0: The regression estimate of the true number of null hypotheses
% - T0sd: The standard deviation associated with our estimate of T0
%
% Schweder T and Spjotvoll E (1982) "Plots of P-values to evaluate many
% tests simultaneously" Biometrika, 69(3):493-502
%
% RTB wrote it, post-150 mile ride on the fixer, 20 May 2012

% For some fun data:
% load prostatecancerexpdata
% pvalues = mattest(dependentData, independentData, 'permute', true);
% [T0] = PvalPlot(pvalues);
% compare with MATLAB's Bioinformatics toolbox:
% [fdr, q] = mafdr(pvalues, 'showplot', true);

if nargin < 2, pFlag = 0; end

p = sort(unique(P),1,'descend');
Np = ones(length(p),1) .* NaN;

for k = 1: length(p)
    Np(k) = sum(P > p(k));
end
figure; plot(1-p,Np,'k.'); hold on;
xlabel('1-p');
ylabel('N_p');
title('P-value plot according to Schweder & Spjotvoll 1982');

% calculate T0 at some reasonable value of p
p0 = 0.5;
T0 = sum(P > p0) / (1-p0);

% I'm still not sure about T0sd; see p. 497 of Schweder & Spjotvoll (1982)
% I get the correct value of 'a', but I don't come up with the same value
% of var(T0) as they do in the text (=8.9, for p0=0.3 and Np = 18).
pn = [0.5 -0.5 -T0];
r = roots(pn);       % solving for 'a' in T0 = 0.5a*(a-1)
a = r(1);
% The value of 0.522 comes from Table 1, p. 497 of S&S 1982 for a
% correlation of 0.5 and 1-p = 0.7
VarNp = (0.5*a*(a-1)*p0*(1-p0)) + (a*(a-1)*(a-2) * (0.522 - ((1-p0)^2)));
VarT0 = VarNp / ((1-p0)^2);
T0sd = sqrt(VarT0);

% Fit a line to the lower part of the curve
if pFlag
    t = find((1-p) < 0.5);    % select the straight part of the curve
    a = polyfit(1-p(t),Np(t)',1);    % a(1) is slope; a(2) is y-intercept
    h = line([0 1],[a(2) a(1)+a(2)]);
    set(h, 'LineStyle', '-', 'LineWidth',2,'Color', 'r');
end
