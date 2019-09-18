% plotMany.m
%
% estimate GVT bias as a function of n and k
% see GVTsim.m for details
% adding silly comment

allN = [50,100,200,500];
allK = [1,2,3];

allDmean = zeros(length(allK),length(allN));
allDstd = zeros(length(allK),length(allN));
allDmedian = zeros(length(allK),length(allN));
for k = 1:length(allK)
    for n = 1:length(allN)
        [D,~,~] = GVTsim(allN(n),allK(k),0.5,10000,0,0);
        allDmean(k,n) = mean(D,'omitnan');
        allDstd(k,n) = std(D,'omitnan');
        allDmedian(k,n) = median(D,'omitnan');
    end
end

figure
hold on
% What I want to do is get the medians plotted in the same color, but as
% dashed lines. This feels like a kluge:
legStr = cell(1,length(allK));
for k = 1:length(allK)
    % error bars are huge
    %hp(k) = errorbar(allN,allDmean(k,:),allDstd(k,:),'o-');
    hp(k) = plot(allN,allDmean(k,:),'o-');
    plot(allN,allDmedian(k,:),'o--','Color',get(hp(k),'Color'));
    legStr{k} = ['k = ',num2str(allK(k))];
end
legend(hp,legStr,'Location','southeast')
xlabel('N');
ylabel('Bias');
title('Counting bias in GVT 1985');
%set(gca,'FontSize',12);