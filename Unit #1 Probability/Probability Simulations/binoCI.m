% binoCI.m: bootstrapping binomial confidence intervals
%
% RTB wrote it, 28 Oct. 2017, reading Nate Silver's book, Signal/Noise

% "During a given season, a true .275 hitter has about a 10 percent chance
% of hitting .300 and a 10 percent chance of hitting .250 on the basis of
% luck alone." p. 79 (footnote 8 says this is based on 500 at-bats in a season)

% What does Nate mean here? Well he's probably talking about the 80%
% confidence interval for the binomial proportion of 138 of 500:
nABs = 500;
trueBA = 0.275;
nHits = round(trueBA*nABs);
myAlpha = 0.2;  % for 80% CI
% Method #1
[pHat,pCI] = binofit(nHits,nABs,myAlpha);
% Yup. pCI = 0.2501    0.3033

% Let's bootstrap this. Key concept is to generate the "raw" data:
x = [ones(nHits,1);zeros(nABs-nHits,1)];    % raw data for season

% anonymous function to calculate a proportion:
prop = @(n) sum(n) ./ length(n);

% Method #2: use 'bootci' to do the heavy lifting:
nBoot = 10000;
ci = bootci(nBoot,{prop,x},'alpha',myAlpha,'type','percentile');
% ci = 0.2500   0.3020

% Method #3: we can do it with 'bootstrp'
allProps = bootstrp(nBoot,prop,x);
allPropsSorted = sort(allProps);
ciLo = allPropsSorted(floor((myAlpha/2 * nBoot)));
ciHi = allPropsSorted(nBoot - floor((myAlpha/2 * nBoot)));
% ciLo = 0.250, ciHi = 0.302

% Method #4: really old-fashioned way, with a 'for' loop and 'unidrnd':
allPropsBS = zeros(nBoot,1);
for k = 1:nBoot
    allPropsBS(k) = prop(x(unidrnd(nABs,nABs,1)));
end
allPropsSorted = sort(allPropsBS);
ciLo = allPropsSorted(floor((myAlpha/2 * nBoot)));
ciHi = allPropsSorted(nBoot - floor((myAlpha/2 * nBoot)));
% ciLo = 0.250, ciHi = 0.302

% Method #5: the parametric bootstrap:
allPropsPBS = binornd(nABs,trueBA,nBoot,1) ./ nABs;
allPropsSorted = sort(allPropsPBS);
ciLo = allPropsSorted(floor((myAlpha/2 * nBoot)));
ciHi = allPropsSorted(nBoot - floor((myAlpha/2 * nBoot)));
% ciLo = 0.250, ciHi = 0.300