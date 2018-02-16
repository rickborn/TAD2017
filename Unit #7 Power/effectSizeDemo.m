% effectSizeDemo.m 

% How likely are we to get a given effect size under H0?
%
% Inspired by a failure-to-listen of Slater Sharp
%
% RTB wrote it, 19 Oct. 2017, gorgeous autumn day

%% Distribution of d' under H0

% Assume H0 is true. Draw 'n' samples from same distribution and measure
% dPrime. How is it distributed? What is the probability of getting a value
% of dPrime > 0.57 by chance?

dPrimeCrit = 0.57;
nSims = 10000;
allN = [5,10,20,50];
pVals = zeros(size(allN));

rng shuffle
figure('Name','Distribution of d-prime under H0');
for k = 1:length(allN)
    subplot(2,2,k);
    allSamples = randn(allN(k),2,nSims);
    allDprime = diff(mean(allSamples));
    histogram(allDprime);
    xlabel('dPrime'); ylabel('#');
    title(['N=' num2str(allN(k)) ', H0 is True']);
    
    if k == 1
        maxAx = axis;
    end
    ax = axis;
    axis([maxAx(1),maxAx(2),ax(3),ax(4)]);
    
    % What is the probability of getting a d-prime greater than a given value?
    pVals(k) = sum(abs(allDprime) >= dPrimeCrit) / nSims;
    text(maxAx(1)+1,(0.9*(ax(4)-ax(3))+ax(3)),...
        ['p(d''>' num2str(dPrimeCrit) ') = ' num2str(pVals(k))]);
end

%% What if there is a real effect, but it is small
dPrimeSim = 0.1;

nSims = 10000;
allN = [5,10,20,50];
pVals = zeros(size(allN));

rng shuffle
figure('Name','Distribution of d-prime with a small effect');
for k = 1:length(allN)
    subplot(2,2,k);
    S1 = normrnd(1,1,allN(k),nSims);
    S2 = normrnd(1+dPrimeSim,1,allN(k),nSims);
    allDprime = mean(S2) - mean(S1);
    histogram(allDprime);
    xlabel('dPrime'); ylabel('#');
    title(['N =' num2str(allN(k)) ', True d'' = ' num2str(dPrimeSim)]);
    
    if k == 1
        maxAx = axis;
    end
    ax = axis;
    axis([maxAx(1),maxAx(2),ax(3),ax(4)]);
    
    % What is the probability of getting a d-prime greater than a given value?
    pVals(k) = sum(abs(allDprime) >= dPrimeCrit) / nSims;
    text(maxAx(1)+1,(0.9*(ax(4)-ax(3))+ax(3)),...
        ['p(d''>' num2str(dPrimeCrit) ') = ' num2str(pVals(k))]);
end
