function [FPrate] = dfSim2(nInit,nAddObs,nMax,myAlpha,nSims,pFlag)

% dfSim2: Simulation of expected false positive rates when data collection
% ends upon obtaining significance, as a function of the frequency with
% which significance tests are performed. This replicates figure 1 from
% Simmons et al. 2011.
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
% RTB wrote it, winter rain storm, 21 Dec 2012; Gill, MA

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
rng('shuffle');

if nargin < 6, pFlag= 0; end
if nargin < 5, nSims = 1000; end
if nargin < 4, myAlpha = 0.05; end
if nargin < 3, nMax = 50; end
if nargin < 2, nAddObs = [1,2,5,10,20]; end
if nargin < 1, nInit = 10; end

FPrate = ones(length(nInit),length(nAddObs)) .* NaN;

for jVal = 1:length(nInit)
    for iVal = 1:length(nAddObs)
        FP = zeros(nSims,1);
        for jSim = 1:nSims
            allSims = randn(nMax,2);
            allNdx = [nInit(jVal):nAddObs(iVal):nMax];
            for k = 1:length(allNdx)
                if ttest2(allSims([1:allNdx(k)],1), allSims([1:allNdx(k)],2), ...
                        'Alpha',myAlpha)
                    FP(jSim) = 1;
                    break
                end
            end
            FPrate(jVal,iVal) = (sum(FP) / nSims) * 100;
        end
    end
end

if pFlag
    figure, plot(nAddObs,FPrate,'b-'); hold on;
    hp = plot(nAddObs,FPrate,'ko');
    set(hp,'MarkerFaceColor','k');
    ax = axis;
    axis([ax(1), ax(2), 0, (floor(ax(4)/5)+1)*5]);
    hl = line([ax(1),ax(2)], [myAlpha*100,myAlpha*100]);
    set(hl,'Color','r','LineStyle','--');
    xlabel('Number of Additional Observations Before Performing Another t-Test');
    ylabel('Percentage of False-Positive Results');
    x = [nAddObs;nAddObs];
    text(x(:)+0.5,FPrate(:)+0.5,num2str(FPrate(:)));
end