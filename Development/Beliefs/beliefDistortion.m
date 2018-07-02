function [subjP] = beliefDistortion(trueP,q)

% q(1) = p0, or the anchor point
% q(2) = gamma, or the slope

Q = (q(2).* logit(trueP)) + ((1 - q(2)).*logit(q(1)));
subjP = 1 ./ (1 + exp(-Q));


end

function [loP] = logit(p)

loP = log(p ./ (1-p));

end