function [pFinal] = RunningProportion2(nRolls)

% RunningProportion2.m: rolls of a die

if nargin < 1, nRolls = 1000; end

nSides = 6;
rollSequence = unidrnd(nSides,nRolls,1);
allN = [1:nRolls]';
runProp = cumsum(rollSequence) ./ allN;
pFinal = runProp(end);
expectedVal = sum((ones(1,nSides) ./ nSides) .* [1:nSides]);

figure
semilogx(allN,runProp,'o-');
xlabel('Roll number');
ylabel('Running average');
ax = axis;
line([ax(1),ax(2)],[expectedVal,expectedVal],'Color','k','LineStyle','--');
text(ax(2)/100,((ax(4)-ax(3))/2)+ax(3),['Mean = ',num2str(pFinal)]);
title('The law of large numbers');
end

