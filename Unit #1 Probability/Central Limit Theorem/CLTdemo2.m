% CLTdemo2.m
%
% This version of the demo explicitly compares the SEM you get when
% calculating by the formula (SDsamp / sqrt(N)) versus the SD of the
% sampling distribution of the mean.
%
% RTB wrote it, 03 September 2020, NEH, Maine

clearall
rng shuffle

R = 5;
nSum = [4,16,25,100];
nSim = 1000;
whichDist = 'unidrnd';

semMean = zeros(length(nSum),1);
sdMean = zeros(length(nSum),1);
muMean = zeros(length(nSum),1);
ctr = 0;
figure('Name','Sampling Distribution of the Mean');
for thisSum = nSum
    ctr = ctr+1;
    allSims = unidrnd(R,[thisSum,nSim]);
    allMeans = mean(allSims);
    
    % formula to calculate the average SEM across columns
    semMean(ctr) = mean(std(allSims) ./ sqrt(thisSum));
    
    % SD of the sampling distribution of the mean, which is, by definition,
    % the standard error of the mean:
    sdMean(ctr) = std(allMeans);
    muMean(ctr) = mean(allMeans);
    
    subplot(2,2,ctr);
    histogram(allMeans);
    hold on
    ax = axis;
    line([muMean(ctr),muMean(ctr)],[ax(3),ax(4)],'Color','r');
    tStr = sprintf('R = %d, SD = %.2f',R,sdMean(ctr)); title(tStr);
    xStr = sprintf('Mean of %d random draws',thisSum); xlabel(xStr);
    ylabel('#');
      
end

figure('Name','SEM Two Ways');
plot(semMean,sdMean,'ro');
xlabel('SEM formula');
ylabel('SD of Sampling Dist of Mean');
title('Comparison of Std. Error of the Mean');
hold on
axis square
ax = axis;
xyMin = min([ax(1),ax(3)]);
xyMax = min([ax(2),ax(4)]);
line([xyMin,xyMax],[xyMin,xyMax]);