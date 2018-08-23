function counts = countSpikes(spikes, timeStepS, startS, endS)

if (nargin < 4)
    endS = length(spikes) * timeStepS;
end
if (nargin < 3)
    startS = 0;
end
trains = size(spikes, 1);
counts = zeros(1, trains);
startBin = startS / timeStepS + 1;
endBin = floor(endS / timeStepS);

for train = 1:trains
    counts(train) = sum(spikes(train, startBin:endBin));
end
end