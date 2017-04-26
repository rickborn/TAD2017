function [rawImage] = mapPinData(pwaNum,N)

% mapPinData.m: Takes raw PIN values and creates 2D heatmap
% Reproduces graphic at http://www.datagenetics.com/blog/september32012/index.html
%
% Inputs:
% - pwaNum: vector of PIN numbers
% - N: range over which to analyze (default = 99)
%
% Outputs:
% - rawImage: log10 of frequency counts for each grid position
%
% RTB wrote it, 25 August 2015 (bored on plane ride, MKE to BOS)

if nargin < 2, N = 99; end
if N < 10, N = 10; end
if N > 99, N = 99; end

% NOTE: If the number of PINs if very large, may want to subsample
% y = datasample(pwaNum,100000,1);

% converting is slow, so better to do this once outside the function
if isstr(pwaNum), pwaNum = str2num(pwaNum); end

% string of pairs of digits from 00 to 99
allNumStrs = ['00';'01';'02';'03';'04';'05';'06';'07';'08';'09';[num2str([10:N]')]];
rawImage = ones(size(allNumStrs,1),size(allNumStrs,1)) .* NaN;
for iX = 1:size(allNumStrs,1)       % column
    for iY = 1:size(allNumStrs,1)   % row
        thisPIN = str2num([allNumStrs(iX,:),allNumStrs(iY,:)]);
        nMatches = sum(pwaNum == thisPIN);
        rawImage(end-iY+1,iX) = log10(nMatches);
    end
end

figure, imagesc(rawImage);
colormap(hot); hCB=colorbar;

% Make it pretty
xlabel('Left two digits','FontSize',18); 
ylabel('Right two digits','FontSize',18);
title('Grid Plot of PIN Data','FontSize',20);
% Because we start with '0000' in the lower left corner, the tick labels
% are off by one and upside down (y-axis)
originalXTickNames = get(gca,'XTickLabel');
newXTickNames = num2str(str2num(char(originalXTickNames)) - 1);
set(gca,'XTickLabel',newXTickNames);
originalYTickNames = get(gca,'YTickLabel');
newYTickNames = num2str(100 - str2num(char(originalYTickNames)));
set(gca,'YTickLabel',newYTickNames);

% make ticks point outwards
set(gca,'TickDir','Out');

% prevent ticks from being re-labeled when we change size
set(gca,'XTickMode','manual');
set(gca,'YTickMode','manual');
axis('image');

% erase tick labels
% set(gca,'XTickLabel','');
% set(gca,'YTickLabel','');

% erase ticks
% set(gca,'XTick',[]);
% set(gca,'YTick',[]);