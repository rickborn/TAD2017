% npPlot.m: plots alpha and beta for different criteria
%
% attempt to reprodcue fig. 2.2 of Efron

allCrit = [1,0.8,0.6,0.4,0.2,0,-0.2,-0.4,-0.6,-0.8,-1];
allAlpha = zeros(size(allCrit));
allBeta = zeros(size(allCrit));
allSamp = [10,20,30];
allH = zeros(size(allSamp));
allLeg = cell(1,3);
cStr = ['rbk'];
dPrime = 0.57;
% nSamp = 10;

figure

for n = 1:length(allSamp)
    for k = 1:length(allCrit)
        [allAlpha(k),allBeta(k)] = npSim(allCrit(k),dPrime,allSamp(n),10000);
    end
    
    allH(n) = plot(allAlpha,allBeta,[cStr(n),'o']);
    allLeg{n} = ['n=',num2str(allSamp(n))];
    hold on
    plot(allAlpha,allBeta,'k-');
    
end

xlabel('\alpha','FontWeight','bold');
ylabel('\beta','FontWeight','bold');
title('Neyman-Pearson Hypothesis Testing');
set(gca,'FontSize',14);
legend(allH,allLeg);
ax = axis;
%text((ax(2)-ax(1))/2,(ax(4)-ax(3))/2,['n=', num2str(nSamp) ', d''=',num2str(dPrime)]);