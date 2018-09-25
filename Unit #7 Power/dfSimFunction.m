function [FPrate] = dfSimFunction(nInit,nAddObs,nMax,myAlpha,nSims)

% dfSimFunction: Simulation of expected false positive rates when data
% collection ends upon obtaining significance, as a function of the
% frequency with which significance tests are performed. This replicates
% figure 1 from Simmons et al. 2011.
%
% Inputs:
% - nInit, # of initial observations prior to first t-test (default, 10)
% - nAddObs, # of additional observations added before next t-test
% - nMax, maximum # of observations to make (default, 50)
% - myAlpha, criterion for rejecting H0 (default, 0.05)
% - nSims, # of simulations to run (default, 1000)
% 
% Outputs:
% - FPrate, percentage of simulations yielding a false positive
%
% RTB wrote it, winter rain storm, 21 Dec 2012; Gill, MA
% RTB converted to a simpler function 19 Sept. 2018 (originally dfSim2.m)
% in order to share with TAD (NB 308qc) class on week #7.

% Reference:
% Simmons JP, Nelson LD, Simonsohn U. False-positive psychology:
% undisclosed flexibility in data collection and analysis allows presenting
% anything as significant. Psychol Sci. 2011 Nov;22(11):1359-66.

% NOTE: I originally wrote SampleSizeDF.m to replicate figure 2 of Simmons
% et al., then built this simulation, which can replicate fig. 1, off of that.

% re-init random number generator
rng shuffle;

if nargin < 5, nSims = 1000; end
if nargin < 4, myAlpha = 0.05; end
if nargin < 3, nMax = 100; end
if nargin < 2, nAddObs = 10; end
if nargin < 1, nInit = 20; end

% The logic for this simulation is that we generate all nMax pairs of
% observations at once. Then, in our inner 'for' loop where we run the
% t-test, we just extend the the length of our index into the rows of
% 'allSims'.
FP = zeros(nSims,1);
allNdx = [nInit:nAddObs:nMax];
for jSim = 1:nSims
    allSims = randn(nMax,2);
    for k = 1:length(allNdx)
        if ttest2(allSims(1:allNdx(k),1), allSims(1:allNdx(k),2),'Alpha',myAlpha)
            FP(jSim) = 1;
            break
        end
    end
end
FPrate = (sum(FP) / nSims) * 100;