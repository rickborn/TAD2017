function [T0,T0sd] = PvalPlot(P,pCrit,pFlag)

% PVALPLOT.M: Plot of P-values according to Schweder & Spjotvoll 1982
%
% [T0] = PvalPlot(P,pCrit,pFlag)
%
% Inputs:
% - P: a columns vector of p-values from a large number of comparisons
% - pCrit: p-value at which to calculate the slope
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

% For some data:
% load prostatecancerexpdata
% pvalues = mattest(dependentData, independentData, 'permute', true);
% [T0,T0sd] = PvalPlot(pvalues,0.3,1);
% compare with MATLAB's Bioinformatics toolbox:
% [fdr, q] = mafdr(pvalues, 'showplot', true);
%
% How should one think about T0? Since there are about T0 true null hypotheses we
% should, when aiming at an overall level a for at least one false
% rejection, use level alpha/T0 for the individual tests. This is an
% improvement over the same Bonferroni argument applied to all tests which
% would have given level alpha/N for the individual tests. In the case of
% the prostate cancer data, we get T0 = 16,981, which is less conservative
% than N = 22,215.

% Also try on Efron's smaller prostate cancer data set:
% load prostateData.mat
% pvalues = mattest(patients, controls, 'permute', true);
% [T0,T0sd] = PvalPlot(pvalues,0.3,1);

if nargin < 3, pFlag = 0; end
if nargin < 2, pCrit = 0.3; end
if size(P,1) == 1, P = P'; end
if pCrit > 0.9, pCrit = 0.3; end

p = sort(unique(P),1,'descend');
Np = ones(length(p),1) .* NaN;

for k = 1: length(p)
    Np(k) = sum(P > p(k));
end
figure; plot(1-p,Np,'k.'); hold on;
xlabel('1-p');
ylabel('N_p (# p-vals > p)');
title('P-value plot according to Schweder & Spjotvoll 1982');

% calculate T0 at some reasonable value of p
T0 = sum(P > pCrit) / (1-pCrit);

% There are different ways of computing the variance in our estimate of T0.
% The first one in S&S (section 3.3, p. 496-7) is not appropriate for array
% data, since it was accounting for larger correlations introduced by doing
% pairwise comparisons between 17 means (i.e. 17-choose-2 = 136 t-tests).
% The more appropriate formula is the one on p. 498.
pn = [0.5 -0.5 -T0];
r = roots(pn);       % solving for 'a' in T0 = 0.5a*(a-1)
a = r(1);
% The value of 0.522 comes from Table 1, p. 497 of S&S 1982 for a
% correlation of 0.5 and 1-p = 0.7
% varNp = (0.5*a*(a-1)*pCrit*(1-pCrit)) + (a*(a-1)*(a-2) * (0.522 - ((1-pCrit)^2)));
varNp = (0.5*a*(a-1)*pCrit*(1-pCrit));
varT0 = varNp / ((1-pCrit)^2);
T0sd = sqrt(varT0);

% Fit a line to the lower part of the curve
if pFlag
    t = find(1-p < (1-pCrit));    % select the straight part of the curve
    a = polyfit(1-p(t),Np(t),1);    % a(1) is slope; a(2) is y-intercept
    h = line([0 1],[a(2) a(1)+a(2)]);
    set(h, 'LineStyle', '-', 'LineWidth',2,'Color', 'r');
    hold on
    % add lines at +/- 2 * T0sd:
    h = line([0 1],[a(2)+(2*T0sd), a(1)+a(2)+(2*T0sd)]);
    set(h, 'LineStyle', '--', 'LineWidth',1,'Color', 'r');
    h = line([0 1],[a(2)-(2*T0sd), a(1)+a(2)-(2*T0sd)]);
    set(h, 'LineStyle', '--', 'LineWidth',1,'Color', 'r');
end
