function [FPrate,nSimsActual] = dfSim3(nInit,nAddObs,nMax,myAlpha,pCrit,nSims,pFlag)

% dfSim3: Simulation of expected false positive rates when data collection
% ends upon obtaining significance, as a function of the frequency with
% which significance tests are performed. This replicates figure 1 from
% Simmons et al. 2011.
%
% In this version, I include another parameter, 'pCrit', to try to make the
% simulation slighly more real. If the first t-test does not look
% "promising" (i.e. p > pCrit), then we don't bother to launch into our
% routine, and we don't count this as a simulation.
%
% Inputs:
% - nInit, # of initial observations prior to first t-test (default, 10)
% - nAddObs, # of additional observations added before next t-test
% - nMax, maximum # of observations to make (default, 50)
% - myAlpha, criterion for rejecting H0 (default, 0.05)
% - pCrit, critical p-value for continuing after the first test (default,0.2)
% - nSims, # of simulations to run (default, 1000)
% - pFlag: if 1, plot results
% 
% Outputs:
% - FPrate, percentage of simulations yielding a false positive
%
% RTB wrote it, winter rain storm, 21 Dec 2012; Gill, MA
% RTB modified to 'dfSim3.m', 30 Sept 2017; rainy Saturday at homoe in JP

% Reference:
% Simmons JP, Nelson LD, Simonsohn U. False-positive psychology:
% undisclosed flexibility in data collection and analysis allows presenting
% anything as significant. Psychol Sci. 2011 Nov;22(11):1359-66.

% NOTE: I originally wrote SampleSizeDF.m to replicate figure 2 of Simmons
% et al., then built this simulation, which can replicate fig. 1, off of that.

% Execultion speed:
% tic;FPrate = dfSim2(10,[1,2,5,10,20],50,0.05,1000,1);toc
% Elapsed time is 37.2 seconds on my Lenovo T61
% Elapsed time is 28.0 seconds on my Dell Latitude E7440

% re-init random number generator
rng shuffle;
%rng default;

if nargin < 7, pFlag= 0; end
if nargin < 6, nSims = 1000; end
if nargin < 5, pCrit = 0.2; end
if nargin < 4, myAlpha = 0.05; end
if nargin < 3, nMax = 50; end
if nargin < 2, nAddObs = [2,5,10,20]; end
if nargin < 1, nInit = 10; end

FPrate = ones(length(nInit),length(nAddObs)) .* NaN;

% The logic for this simulation is that we generate all nMax pairs of
% observations at once. Then, in our innter 'for' loop where we run the
% t-tests, we just extend the the length of our index into 'allSims'.
for jVal = 1:length(nInit)
    for iVal = 1:length(nAddObs)
        FP = zeros(nSims,1);
        nSimsActual = 0;
        for jSim = 1:nSims
            allSims = randn(nMax,2);
            allNdx = [nInit(jVal):nAddObs(iVal):nMax];
            [~,p] = ttest2(allSims([1:allNdx(1)],1), allSims([1:allNdx(1)],2), ...
                        'Alpha',myAlpha);
            if p < pCrit
                nSimsActual = nSimsActual + 1;
            for k = 2:length(allNdx)
                if ttest2(allSims([1:allNdx(k)],1), allSims([1:allNdx(k)],2), ...
                        'Alpha',myAlpha)
                    FP(jSim) = 1;
                    break
                end
            end
            end
            FPrate(jVal,iVal) = (sum(FP) / nSimsActual) * 100;
        end
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