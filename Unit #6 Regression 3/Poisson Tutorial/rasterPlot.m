function rasterPlot(spikes, timeStepS)

figure(1);

times = [0:timeStepS:timeStepS * (length(spikes) - 1)];
axes('position', [0.1, 0.1, 0.8, 0.8]);
axis([0, length(spikes) - 1, 0, 1]);
trains = size(spikes, 1); 
ticMargin = 0.01;                                      % gap between spike trains (full scale is 1)
ticHeight = (1.0 - (trains + 1) * ticMargin) / trains;

for train = 1:trains
    spikeTimes = find(spikes(train, :) == 1);
    yOffset = ticMargin + (train - 1) * (ticMargin + ticHeight);
    for i = 1:length(spikeTimes)
        line([spikeTimes(i), spikeTimes(i)], [yOffset, yOffset + ticHeight]);
    end
end

xlabel('Time (msec)')
title('Raster plot of spikes');
end