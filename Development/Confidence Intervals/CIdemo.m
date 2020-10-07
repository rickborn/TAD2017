function [trueCIinex,mitCIinex] = CIdemo(nSamp,alpha,nSims,pFlag)

% CIdemo.m: Compares two definitions of a frequentist confidence interval
%
% e.g. [trueCIinex,mitCIinex] = CIdemo(10,0.2,1000,1);
%
% Inputs:
% - nSamp: # of samples in each experiment (default = 10)
% - alpha: defines (1 - alpha)*100 % CI (default = 0.2 for 80% CI)
% - nSims: # of simulations to run (default = 10,000)
% - pFlag: 1 = plot results (default = 0 for no plot)
%
% Outputs:
% - trueCI: # of sims on which true mean is included/excluded from CI
% - mitCI: # of sims on which 1st CI includes/doesn't include repeats of the mean
%
% RTB wrote it, 29 September 2020, NE Harbor, grey, foggy morning

% cd 'C:\usr\rick\doc\Committees\PIN\PIN Director\Courses\Stats\TAD\TAD Code\Development\Confidence Intervals'

% This demo was inspired by a mistake in a beautiful study at MIT on the
% effects of college scholarships. See:
% https://evaluatingcollegesupport.mit.edu/
%
% What's cool about the study is that they used a randomized design to test
% for causality. But in their legends, their explanation for the 95% CI
% caught my eye:
%
% "The small interval markers overlapping the top of each bar indicate a 95
% percent confidence interval. This means that if the STBF scholarship
% study were repeated 100 times with different random samples of students,
% 95 of those 100 studies would produce results within the interval area."
%
% The final clause is where they fall down. The 95% CI is itself a random
% variable: for each one of those 100 repeats, the calculated 95% CI would
% be different. We would expect that about 95% of those CIs would contain
% the true mean. There is nothing special about the 1 CI we happened to
% calculate. BUT, the erroneous statement--though wrong in principle--might
% not actually be so far off in practice. Let's check; and then let's see
% how it is affected by the sample size, 'nSamp'.

% CONCLUSION: The MIT definition is not only wrong in principle, but very
% wrong in practice. It will occasionally get it right, but it is often
% disastrously wrong, even for larg N.

% We will ask two questions:
% 1. How often does our calculated CI contain the true mean? (proper CI)
% 2. How often does our "first" CI contain subsequent means? (MIT def'n)

% defaults:
if nargin < 4, pFlag = 0; end
if nargin < 3, nSims = 10000; end
if nargin < 2, alpha = 0.2; end
if nargin < 1, nSamp = 10; end

% pre-allocate our arrays:
allMeans = zeros(nSims,1);
allCI = zeros(2,nSims);

% So, let's pretend our first CI is "special"; i.e. that it is the gold
% standard to test against the means of future repeats. For this, and for
% our subsequent simulated repeats, we will just draw 'nSamp' samples from
% a standard normal distribution and use 'normfit' to calculate our means
% and confidence intervals.
x = randn(nSamp,1);
[~,~,muCI,~] = normfit(x,alpha);

for k = 1:nSims
    [allMeans(k),~,allCI(:,k),~] = normfit(randn(nSamp,1),alpha);
end

% 1. How often does our calculated CI contain the true mean? (Def. #2)
nTrueInclude = sum(allCI(1,:) < 0 & allCI(2,:) > 0);
nTrueExclude = nSims - nTrueInclude;
trueCIinex = [nTrueInclude,nTrueExclude];

% 2. How often does our "first" CI contain subsequent means? (Def. #1, MIT)
nMITinclude = sum(allMeans > muCI(1,1) & allMeans < muCI(2,1));
nMITexclude = nSims - nMITinclude;
mitCIinex = [nMITinclude,nMITexclude];

% plot side-by-side bar graphs of the results:
if pFlag
    subplot(1,2,2);
    X = categorical({'Contains mu','Excludes mu'});
    X = reordercats(X,{'Contains mu','Excludes mu'});
    bar(X,[nTrueInclude,nTrueExclude]);
    hold on
    % put a marker at where each SHOULD be:
    yShould = [nSims * (1 - alpha), nSims * alpha];
    hp = plot(X,yShould,'rs');
    set(hp,'MarkerFaceColor','r');
    ylabel('#');
    %title(['Proper def. of ' num2str((1-alpha)*100,2) '% CI']);
    title(['Explanation #2 of ' num2str((1-alpha)*100,2) '% CI']);
    ylim([0,nSims]);
    hold off
    
    subplot(1,2,1);
    bar(X,[nMITinclude,nMITexclude]);
    hold on
    hp = plot(X,yShould,'rs');
    set(hp,'MarkerFaceColor','r');
    ylabel('#');
    % title(['MIT def. of ' num2str((1-alpha)*100,2) '% CI']);
    title(['Explanation #1 of ' num2str((1-alpha)*100,2) '% CI']);
    ylim([0,nSims]);
    hold off
end