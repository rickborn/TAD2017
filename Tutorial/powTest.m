% powTest.m

nA = 25;
nB = 25;
effSize = 0.5;

nSim = 100000;
allSig = zeros(nSim,1);

for k=1:nSim
    a = normrnd(0,1,nA,1);
    b = normrnd(effSize,1,nB,1);
    allSig(k) = ttest2(a,b);
end

% find false negatives:
pFN = (length(allSig) - (sum(allSig))) / nSim;
powEst = 1 - pFN;

% calculate power by the formula:
powTrue = sampsizepwr('t2',[0,1],0.5,[],nA);

