% pCurveFamily.m
%
% Generates range of p-curves
% Reproduces Figure 1 of Simonsohn et al. 2014
%
% RTB wrote it, some time in 2015

allN = [20,50,100,20,50,100,20,50,100];
allD = [0.42,0.26,0.18,0.64,0.4,0.28,0.91,0.57,0.4];
nSim = 100000;

for k = 1:9
    subplot(3,3,k);
    [~,~,~] = pCurve(allD(k),allN(k),nSim);
end