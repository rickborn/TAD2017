function [pFinal] = RunningProportion(nFlips,pHeads)

% RunningProportion.m: sequence of coin flips

if nargin < 2, pHeads = 0.5; end
if nargin < 1, nFlips = 500; end

flipSequence = rand(nFlips,1) < pHeads;
allN = [1:nFlips]';
runProp = cumsum(flipSequence) ./ allN;
pFinal = runProp(end);

figure
semilogx(allN,runProp,'o-');
xlabel('Flip number');
ylabel('Proportion heads');
ax = axis;
line([ax(1),ax(2)],[pHeads,pHeads],'Color','k','LineStyle','--');
text(ax(2)/100,((ax(4)-ax(3))/2)+ax(3),['End Proportion = ',num2str(pFinal)]);

end

