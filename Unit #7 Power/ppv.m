function [y] = ppv(R,alpha,power)

% ppv.m: positive predictive value
%
% Inputs:
% - R: prior odds in favor of a non-null result (default = 0.25)
% - alpha: type I error (default = 0.05)
% - power: 1 - type II error (default = 0.8)

if nargin < 3, power = 0.8; end
if nargin < 2, alpha = 0.05; end
if nargin < 1, R = 0.25; end

y = (power.*R) ./ (power.*R + alpha);

end

