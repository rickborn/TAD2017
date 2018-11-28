function [y] = ppv(R,alpha,power)

% ppv.m: positive predictive value
%
% PPV = ppv(R,alpha,power);
%
% Inputs:
% - R: prior odds in favor of a non-null result (default = 0.25)
% - alpha: type I error (default = 0.05)
% - power: 1 - type II error (default = 0.8)
%
% Button KS, Ioannidis JP, Mokrysz C, Nosek BA, Flint J, Robinson ES,
% Munafò MR. Power failure: why small sample size undermines the
% reliability of neuroscience. Nat Rev Neurosci. 2013 May;14(5):365-76.
% doi: 10.1038/nrn3475. Epub 2013 Apr 10. Review. Erratum in: Nat Rev
% Neurosci. 2013 Jun;14(6):451. PubMed PMID: 23571845.
%
% RTB wrote it 17 April 2018 for a slide in his Reproducibility lecture to postdocs

if nargin < 3, power = 0.8; end
if nargin < 2, alpha = 0.05; end
if nargin < 1, R = 0.25; end

y = (power.*R) ./ (power.*R + alpha);

end

