function [FPrate] = dfSim4(nInit,nAddObs,nMax,myAlpha,nSims,pFlag)

% dfSim2: Simulation of expected false positive rates when data collection
% ends upon obtaining significance, as a function of the frequency with
% which significance tests are performed. This replicates figure 1 from
% Simmons et al. 2011. Version 4 was based on dfSim2, but it allows for
% more flexibility in modeling different H0s, such as comparing
% distributions that have identical means but different variances. It will
% be identical to dfSim2 if 'diffVarFlag' is set to 0.
%
% Inputs:
% - nInit, # of initial observations prior to first t-test (default, 10)
% - nAddObs, # of additional observations added before next t-test
% - nMax, maximum # of observations to make (default, 50)
% - myAlpha, criterion for rejecting H0 (default, 0.05)
% - nSims, # of simulations to run (default, 1000)
% - pFlag: if 1, plot results
% 
% Outputs:
% - FPrate, percentage of simulations yielding a false positive
%
% RTB wrote version 4, 12 Feb 2021, in NEH, Maine during the pandemic
% RTB wrote it, winter rain storm, 21 Dec 2012; Gill, MA

% Reference:
% Simmons JP, Nelson LD, Simonsohn U. False-positive psychology:
% undisclosed flexibility in data collection and analysis allows presenting
% anything as significant. Psychol Sci. 2011 Nov;22(11):1359-66.

% NOTE: I originally wrote SampleSizeDF.m to replicate figure 2 of Simmons
% et al., then built this simulation, which can replicate fig. 1, off of that.

% Execultion speed:
% tic;FPrate = dfSim4(10,[1,5,10,20],50,0.05,1000,1);toc
% Elapsed time is 37.2 seconds on my Lenovo T61
% Elapsed time is 28.0 seconds on my Dell Latitude E7440
% Elapsed time is 21.9 seconds on my Dell Ultrabook
% Elapsed time is 9.1 seconds on my MS Surface (15 Dec. 2018)

% re-init random number generator
rng shuffle;
%rng default;

if nargin < 6, pFlag= 0; end
if nargin < 5, nSims = 1000; end
if nargin < 4, myAlpha = 0.05; end
if nargin < 3, nMax = 50; end
if nargin < 2, nAddObs = [1,5,10,20]; end
if nargin < 1, nInit = 10; end

% This is sort of a kluge to allow our two simulated samples to have
% identical means, but different variances. If 'diffVarFlag' is set, the
% second sample will by multiplied by 'varMult'
diffVarFlag = 1;
varMult = 2;
% Increasing the variance of one population does inflate the FP rate, but
% not as much as I expected. The highest rate goes from 22% to 24% with a
% variance multiplier of 2.

% Initialize the variable that will hold our results:
FPrate = ones(length(nInit),length(nAddObs)) .* NaN;

% The logic for this simulation is that we generate all nMax pairs of
% observations at once. Then, in our inner 'for' loop where we run the
% t-tests, we just extend the the length of our index into 'allSims'.
for jVal = 1:length(nInit)
    for iVal = 1:length(nAddObs)
        FP = zeros(nSims,1);
        for jSim = 1:nSims
            % Random draws for both populations:
            allSims = randn(nMax,2);
            % Allow for different variances:
            if diffVarFlag
                allSims(:,2) = allSims(:,2) .* varMult;
            end
            allNdx = [nInit(jVal):nAddObs(iVal):nMax];
            for k = 1:length(allNdx)
                if ttest2(allSims([1:allNdx(k)],1), allSims([1:allNdx(k)],2), ...
                        'Alpha',myAlpha)
                    FP(jSim) = 1;
                    break
                end
            end
        end
        FPrate(jVal,iVal) = (sum(FP) / nSims) * 100;
    end
end

if pFlag
    figure
    cStrs = ['b','r'];
    for k = 1:length(nInit)
        plot(nAddObs,FPrate(k,:),[cStrs(k),'-']); hold on;
        hp = plot(nAddObs,FPrate(k,:),'ko');
        set(hp,'MarkerFaceColor','k');
        
        %x = [nAddObs;nAddObs];
        text(nAddObs(:)-0.5,(FPrate(k,:)+1)',num2str(FPrate(k,:)'));
    end
    ax = axis;
    axis([ax(1), ax(2), 0, (floor(ax(4)/5)+1)*5]);
    hl = line([ax(1),ax(2)], [myAlpha*100,myAlpha*100]);
    set(hl,'Color','r','LineStyle','--');
    xlabel('# of additional observations before performing next test');
    ylabel('False-positive results (%)');
end