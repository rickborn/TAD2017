% tailSim.m: Reich exercise on tail differences
%
% differences in representation in the tails caused by shifts in the mean
% or the variance
%
% RTB wrote it, 17 March 2019; bad head cold, just finished Stephen
% Markley's amazing first novel, "Ohio."

% This was inspired by an exercise in Chapter 11 of David Reich's book on
% ancient DNA ("Who We Are and How We Got Here: Ancient DNA and the New
% Science of the Human Past"). He gives an example of west Africans and
% sprinting: "All the male finalists in the Olympic hundred-meter race
% since 1980, even those from Europe and the Americas, had recent West
% African ancestry." (citing David Epstein's "The Sports Gene: Inside the
% Science of Extraordinary Athletic Performance," New York: Current,
% 2013). He goes on to give two types of explanations for this difference
% at the extremes:

%% Explanation #1: shift in the means

% If, on average, the distribution of sprinting speeds for West Africans
% was shifted by 0.8 SD w/r/t other populations, there would be a 100-fold
% enrichment in the 99.9999999th percentile (which is about 6 SD from the
% mean).
normcdf(6,0.8,1,'upper') / normcdf(6,0,1,'upper')
% ans = 100.9989

% Note that an effect size of 0.8 SD is not a small effect. Cohen (1988)
% would classify it as "large"
% (https://en.wikipedia.org/wiki/Effect_size#Cohen's_d)

%% Explanation #2: difference in variance

% There is an approximately 33% higher genetic diversity in West Africans
% than in Europeans:
normcdf(6,0,sqrt(1.33),'upper') / normcdf(6,0,1,'upper')
% ans = 99.5635

%% So let's look at this over a range of variances

% Some useful numbers:
usPop = 327.2e6;    % Population of the U.S. is 327.2 million
worldPop = 7.7e9;   % Population of the world is about 7.7 billion

% Calibration #'s for height
% (https://tall.life/height-percentile-calculator-age-country/)
mMu = 69.3;     % average US male height in inches
mSD = 2.92;     % standard deviation
fMu = 63.8;
fSD = 2.80;
% effect size
(mMu - fMu) / sqrt(mSD.^2 + fSD.^2)
% ans = 1.36

tailSDs = [3:7];
% 3.1 SD = 99.9th percentile
% 4.2 SD = 99.999th percentile
% 5.2 SD = 99.99999th percentile
% 6.0 SD = 99.9999999th percentile
% As a point of comparison, we would expect there to be only 7.6 people in
% the entire world greater than 6 SD from the mean:
% normcdf(6,0,1,'upper') * worldPop

% About 2800 people in the world are 7 feet tall or taller. Considering
% that the world population is approximately 7.4 billion people, this means
% that the percentage of 7 footers is 0.000038%.
% So, 7 feet is about 5 SD above the mean. Based on the numbers above, we
% can calculate the expected number of people in the world of 7 feet tall
% or taller:
zHgtM = ((7*12) - mMu) / mSD;
nMenPred = normcdf(zHgtM,0,1,'upper') * usPop/2;
zHgtF = ((7*12) - fMu) / fSD;
nWomenPred = normcdf(zHgtF,0,1,'upper') * usPop/2;
nPred = nMenPred + nWomenPred;
% This (i.e. 924) is much smaller than what was reported (2800).
% But note that my height figures were for the U.S., not the world, and I
% wasn't able to find a value for # of people in the U.S. over 7 ft.

% Calibration for differences in variance:
% 2010 meta-analysis of 242 studies found that males have an 8% greater
% variance in mathematical abilities than females.
% 
% Lindberg, Sara M.,Hyde, Janet Shibley,Petersen, Jennifer L.,Linn, Marcia
% C. "New trends in gender and mathematics performance: A meta-analysis."
% Psychological Bulletin, Vol 136(6), Nov 2010, 1123-1135.
% http://psycnet.apa.org/doiLanding?doi=10.1037%2Fa0021276
allVar = [1.08,1.15,1.33,1.5];
allVarLabels = {'1.08','1.15','1.33','1.50'};

% So, for a population with the same mean, but a larger variance, what
% would the relative enrichment be in various extremes of the distribution?
allRatios = zeros(length(allVar),length(tailSDs));

for j = 1:length(allVar)
    for k = 1:length(tailSDs)
        allRatios(j,k) = normcdf(tailSDs(k),0,sqrt(allVar(j)),'upper') / normcdf(tailSDs(k),0,1,'upper');
    end
end

figure
subplot(2,1,1);
semilogy(tailSDs,allRatios,'o-','LineWidth',2);
xlabel('tail above (SD)');
ylabel('fold enrichment');
title('Variance');
legend(allVarLabels,'Location','northwest');

set(gcf,'Position',[456 29 837 779],'Name','Enrichment of the Tails');

%% Same thing, but for different effect sizes
allES = [0.2,0.5,0.8,1.2];
allESlabels = {'0.2','0.5','0.8','1.2'};

% https://en.wikipedia.org/wiki/Effect_size#Cohen's_d
% Very small	0.01	Sawilowsky, 2009
% Small         0.20	Cohen, 1988
% Medium        0.50	Cohen, 1988
% Large         0.80	Cohen, 1988
% Very large	1.20	Sawilowsky, 2009
% Huge          2.0     Sawilowsky, 2009

allRatios = zeros(length(allES),length(tailSDs));

for j = 1:length(allES)
    for k = 1:length(tailSDs)
        allRatios(j,k) = normcdf(tailSDs(k),allES(j),1,'upper') / normcdf(tailSDs(k),0,1,'upper');
    end
end

%figure
subplot(2,1,2);
semilogy(tailSDs,allRatios,'o-','LineWidth',2);
xlabel('tail above (SD)');
ylabel('fold enrichment')
title('Effect Size');
legend(allESlabels,'Location','northwest');

