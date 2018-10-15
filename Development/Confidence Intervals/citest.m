function [p,CI,pwr] = citest(n,dPrime,myAlpha,pFlag)

% citest.m: explores relationship between confidence intervals and signifcance testing
%
% motivated by an article that Brian Healy sent me:
%
% Greenland S, Senn SJ, Rothman KJ, Carlin JB, Poole C, Goodman SN, Altman
% DG. Statistical tests, P values, confidence intervals, and power: a guide
% to misinterpretations. Eur J Epidemiol. 2016 Apr;31(4):337-50. doi:
% 10.1007/s10654-016-0149-3.
%
% RTB wrote it 13 Oct. 2018, cold, raniy autumn day

% Generate our samples:
allSamples = randn(n,2);
% To simulate a real d-prime, we just add our dPrime to one of the samples.
allSamples(:,2) = allSamples(:,2) + dPrime;

CI = zeros(2,2);
% Use one sample t-test to get 95% CI on each population:
[~,~,CI(1,:)] = ttest(allSamples(:,1),0,'Alpha',myAlpha);
[~,~,CI(2,:)] = ttest(allSamples(:,2),0,'Alpha',myAlpha);

% Now run our 2-sample t-test:
[~,p] = ttest2(allSamples(:,2),allSamples(:,1),'Alpha',myAlpha);

% Calculate the power
pwr = sampsizepwr('t2',[0,1],dPrime,[],n);

% Plot stuff:
if pFlag
    %figure
    %boxplot(allSamples);
    hold on
    
    % calculate means and error bars
    % NOTE: The 'errorbar' function plots error bars that are L(i) + U(i)
    % long. That is, it doesn't treat our CI as an interval, but rather as
    % a distance from the mean to the end of each error bar. So we need to
    % subtract each from the mean:
    bothMeans = mean(allSamples);
    ebLo(1,:) = bothMeans' - CI(:,1); % lower error bar
    ebHi(1,:) = CI(:,2) - bothMeans'; % upper error bar
    
    % plot together?
    %he = errorbar([1,2],bothMeans,ebLo,ebHi,'o');
    
    % plot each seperately:
    he1 = errorbar([1],bothMeans(1),ebLo(1),ebHi(1),'s','CapSize',20);
    myBlue = get(he1,'Color');
    set(he1,'LineWidth',2,'MarkerSize',10);
    he2 = errorbar([2],bothMeans(2),ebLo(2),ebHi(2),'s','CapSize',20);
    myRed = get(he2,'Color');
    set(he2,'LineWidth',2,'MarkerSize',10);
    
    % include some relevant text values: n,p-val,power
    %tStr = sprintf('N=%d, dPrime=%.2f, P=%.3f, Pwr=%.2f',n,dPrime,p,pwr);
    tStr = sprintf('P=%.3f',p);
    if p <= myAlpha
        title(tStr,'Color','r');
    else
        title(tStr,'Color','k');
    end
    %text(1.25,2,tStr);
    %xlabel('Group');
    
    % draw some horizontal lines to make comparisons easier:
    axis square
    ax = axis;
    axis([0.75,2.25,ax(3),ax(4)]);
    ax = axis;
    % line at 0 for null
    line([ax(1),ax(2)], [0,0],'Color',myBlue,'LineStyle','--');
    % line at d-prime
    line([ax(1),ax(2)], [dPrime,dPrime],'Color',myRed,'LineStyle','--');
    
    % CI for sample #1:
    %line([1,2],[CI(1,1),CI(1,1)],'Color',myBlue);
    line([1,2],[CI(1,2),CI(1,2)],'Color',myBlue);
    % CI for sample #2:
    line([1,2],[CI(2,1),CI(2,1)],'Color',myRed);
    %line([1,2],[CI(2,2),CI(2,2)],'Color',myRed);
        
%     lStr = sprintf('%d%% CI',100*(1-myAlpha));
%     legend(he,lStr,'Location','North');
    
    set(gca,'XTickLabels',[])
    set(gca,'YTickLabels',[])
    
end
    