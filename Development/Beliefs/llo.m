function [subjP] = llo(trueP,q)

% llo.m: linear log-odds function after Zhang & Maloney 2012
%
% Inputs:
% - trueP: the true underlying frequency or probability
% - q(1) = p0, or the anchor point
% - q(2) = gamma, or the slope
%
% Outputs:
% -subjP: the subjective probability
%
% RTB wrote it, 05 April 2018, plane ride to SF for Jeff's Shukke Tokudo

Q = (q(2).* logit(trueP)) + ((1 - q(2)).*logit(q(1)));
subjP = 1 ./ (1 + exp(-Q));


end

function [loP] = logit(p)

loP = log(p ./ (1-p));

end