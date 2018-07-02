% effectSizeDemo.m 

% How likely are we to get a given effect size under H0?
%
% Inspired by a failure-to-listen of Slater Sharp
%
% RTB wrote it, 19 Oct. 2017, gorgeous autumn day

%% Distribution of d' under H0

% Assume H0 is true. Draw 'n' samples from same distribution and measure
% dPrime. How is it distributed? What is the probability of getting a value
% of dPrime > 0.57 by chance?

dPrimeCrit = 0.57;
nSims = 10000;
allN = [5,10,20,50];
pGreaterThanCrit = zeros(size(allN));
powerSim = zeros(size(allN));
yCrit = -log10(0.05);

rng shuffle
figure('Name','Distribution of d-prime under H0');
for k = 1:length(allN)
    subplot(2,2,k);
    allSamples = randn(allN(k),2,nSims);
    allDprime = squeeze(diff(mean(allSamples)));
    
    % histogram of d-prime values
    yyaxis left
    histogram(allDprime);
    hold on
    xlabel('d-Prime'); ylabel('#');
    title(['N=' num2str(allN(k))]);
    %set(gca,'FontSize',12);
    
    % plot of p-values for a 2-sample t-test
    yyaxis right
    [h,p] = ttest2(squeeze(allSamples(:,1,:)),squeeze(allSamples(:,2,:)));
    h = logical(h);
    plot(allDprime(~h),-log10(p(~h)),'k.');
    plot(allDprime(h),-log10(p(h)),'r.');
    ylabel('-log10(p-Value)');
        
    if k == 1
        %maxAx = axis;
        maxAx = [-2.5,2.5,0,750];
    end
    ax = axis;
    axis([maxAx(1),maxAx(2),ax(3),ax(4)]);
    line([maxAx(1),maxAx(2)],[yCrit,yCrit],'Color','r');
        
    % What is the probability of getting a d-prime greater than a given value?
    pGreaterThanCrit(k) = sum(abs(allDprime) >= dPrimeCrit) / nSims;
    %text(maxAx(1)+1,(0.9*(ax(4)-ax(3))+ax(3)),...
        %['p(d''>' num2str(dPrimeCrit) ') = ' num2str(pGreaterThanCrit(k),2)]);
        
    % What is the probability of getting a false positive?
    powerSim(k) = sum(h) / nSims;
%     text(maxAx(1)+1,(0.9*(ax(4)-ax(3))+ax(3)),...
%         ['p(FP) = ' num2str(pFP(k),2)]);
    
    % What will be the median effect size published?
    pubESright = median(allDprime(h==1 & allDprime' > 0));
    pubESleft = median(allDprime(h==1 & allDprime' < 0));
%     line([pubESright,pubESright],[ax(3),ax(4)],'Color','b','LineStyle','--');
%     line([pubESleft,pubESleft],[ax(3),ax(4)],'Color','b','LineStyle','--');
    plot(pubESright,yCrit,'rs','MarkerFaceColor','r');
    plot(pubESleft,yCrit,'rs','MarkerFaceColor','r');
    line([pubESright,pubESright],[0,yCrit],'Color','r');
    line([pubESleft,pubESleft],[0,yCrit],'Color','r');
end

%% What if there is a real effect?

% Ans. It will still be inflated in the published literature if power is
% low.

% Cohen suggested that d=0.2 be considered a 'small' effect size, 0.5
% represents a 'medium' effect size and 0.8 a 'large' effect size.
dPrimeSim = 0.57;

% Make plots?
pFlag = 1;

if pFlag
    allN = [10,20,50,100];
    s = sprintf('Distribution of d'' with effect size of %0.2f',dPrimeSim);
    figure('Name',s);
else
    allN = [5,10,20,50,100,200,500];
end

nSims = 10000;
pGreaterThanCrit = zeros(size(allN));
powerSim = zeros(size(allN));
pubESright = zeros(size(allN));
yCrit = -log10(0.05);

rng shuffle
for k = 1:length(allN)
    
    allSamples = randn(allN(k),2,nSims);
    % simulate a real d-prime, we just add our mu-offset to one of the
    % samples
    allSamples(:,2,:) = allSamples(:,2,:) + dPrimeSim;
    allDprime = squeeze(diff(mean(allSamples)));    % SD = 1
    
    % Do a one-sided t-test:
    [h,p] = ttest2(squeeze(allSamples(:,2,:)),squeeze(allSamples(:,1,:)),'Tail','Right');
    h = logical(h);
    
    % What is the power?
    powerSim(k) = sum(h) / nSims;
    
    % What will be the median effect size published?
    pubESright(k) = median(allDprime(h==1 & allDprime' > 0));
    
    if pFlag
        subplot(2,2,k);
        % histogram of d-prime values
        yyaxis left
        histogram(allDprime);
        hold on
        xlabel('d-Prime'); ylabel('#');
        s = sprintf('N=%d, True effect size=%0.2f',allN(k),dPrimeSim);
        title(s)
        %title(['N=' num2str(allN(k))]);
        set(gca,'FontSize',12);
        
        % plot of p-values for a 2-sample t-test
        yyaxis right
        plot(allDprime(~h),-log10(p(~h)),'k.');
        plot(allDprime(h),-log10(p(h)),'r.');
        ylabel('-log10(p-Value)');
        
        if k == 1
            %maxAx = axis;
            maxAx = [-2.5,2.5,0,750];
        end
        ax = axis;
        axis([maxAx(1),maxAx(2),ax(3),ax(4)]);
        line([maxAx(1),maxAx(2)],[yCrit,yCrit],'Color','r');
        
        
        text(maxAx(1)+0.2,(0.9*(ax(4)-ax(3))+ax(3)),...
            ['Power = ' num2str(powerSim(k),2)],'FontSize',12);
        
        text(maxAx(1)+0.2,(0.8*(ax(4)-ax(3))+ax(3)),...
            ['Median Published Effect Size = ' num2str(pubESright(k),2)],'FontSize',12);
        
        %pubESleft = median(allDprime(h==1 & allDprime' < 0));
        %     line([pubESright,pubESright],[ax(3),ax(4)],'Color','b','LineStyle','--');
        %     line([pubESleft,pubESleft],[ax(3),ax(4)],'Color','b','LineStyle','--');
        plot(pubESright(k),yCrit,'rs','MarkerFaceColor','r');
        %plot(pubESleft,yCrit,'rs','MarkerFaceColor','r');
        line([pubESright(k),pubESright(k)],[0,yCrit],'Color','r');
        %line([pubESleft,pubESleft],[0,yCrit],'Color','r');
    end
end

relBiasFlag = 0;

% Make a plot of "effect size inflation" as a function of power
figure('Name','Effect Size Inflation');
hold on

if relBiasFlag
    plot(powerSim,(pubESright - dPrimeSim) ./ dPrimeSim,'b*-');
    ylabel('Relative bias');
    ax = axis;
    line([ax(1),ax(2)],[0,0],'Color','r','LineStyle','--');
else
    plot(powerSim,pubESright ./ dPrimeSim,'b*-');
    ylabel('Effect size inflation');
    ax = axis;
    line([ax(1),ax(2)],[1,1],'Color','r','LineStyle','--');
end
xlabel('Power');

%


% Hah! I independently discovered a cool thing. See:
% https://www.nature.com/articles/nrn3475/figures/5
% They express their y-axis as "relative bias", which, in my exaample would
% be: (pubESright - dPrimeSim) ./ dPrimeSim
