%% Multiple Comparisons Exercise
%
% In the world of high throughput molecular biology, we are frequently
% performing many statistical tests on a single data set. In this exercise,
% we will explore the problems inherent in this practice, and introduce
% approaches for solving the problems.
%
% see also: ttest_perm, ttest_norm, DataRand, PvalPlot

% RTB wrote it, 26 May 2012
% YS (Yonatan Savir) adapted for class, 30 May 2012
% RTB made corrections and added 'Step 6' on 19 August 2012

%% Step 1. Read in a large set of Affymetrix data from a prostate cancer study.
% 
% Best CJM, Gillespie JW, Yi Y, Chandramouli GVR, Perlmutter MA, Gathright
% Y, Erickson HS, Georgevich L, Tangrea MA,  Duray PH, Gonzalez S, Velasco
% A, Linehan WM, Matusik RJ, Price DK, Figg WD, Emmert-Buck MR and Chuaqui
% RF (2005) Molecular alterations in primary prostate cancer after androgen
% ablation therapy. Clinical Cancer Research 11:6823–6834.
% 
% The two variables in the MAT-file, dependentData and independentData, are
% two matrices of gene expression values from two experimental conditions,
% where each row corresponds to a gene and each column corresponds to a
% replicate.

load prostatecancerarraydata

%% Step 2. Write a function that will perform a 2-sample t-test on each gene and return a variable containing all of the P-values.
% Create an empty array - allocate space for pValues.
pValues = ones(length(independentData),1) .* NaN;

% A loop over all genes (note that number of genes is larger than the
% number of replicates. H is binary (1-reject 0 can't reject)
for k = 1:length(independentData)
    % TWO sample t-test. Each step compares one row, which is one gene
    [H,pValues(k)] = ttest2(independentData(k,:), dependentData(k,:));
end

% Can this be done without the 'for' loop? Try it!
%[H,pValues] = ttest2(independentData', dependentData');

% How many genes show differential expression?
% (i.e. How many pValues are smaller than 0.05?)
nSig = sum(pValues < 0.05);

% Is something fishy here? (see DataRand.m)
%% Step 3. What kind of distribution do your P-values follow?

figure, hist(pValues);
title('Distribution of p-values from the t-test')
xlabel('p-value')

% What does the corresponding cumulative distribution look like?
x = 0.05:0.1:0.95;          % define our own bin centers
n = hist(pValues,x);        % let 'hist' do the hard part
C = cumsum(n);             
C = cumsum(n)/max(cumsum(n));% set the maximum to 1
figure, plot(x,C);

title('Cumulative distribution of p-values from the t-test')
xlabel('p-value')
ylabel('Percentage of probes')

% What is the slope of this line?
%% Step 4. P-value Plots
% Schweder T and Spjotvoll E (1982) "Plots of P-values to evaluate many
% tests simultaneously" Biometrika, 69(3):493-502
%
% The key to their approach is that, if HO is true, then p-values are
% uniformly distributed over the interval (0,1). Schweder & Spjotvoll used
% these facts to devise a procedure for plotting P-values from a large
% number of simultaneous comparisons. For each P-value, p, obtained, one
% calculates the number of P-values in the set that are greater than that
% P-value, Np, and then plots Np vs. 1-p. From step 6 above you'll
% recognize that Np is just "area beyond" in the cumulative density of the
% P-value distribution. Using this convention, and plotting it versus 1-p
% (instead of p) simply puts the interesting, extreme values (low p) to the
% far right of the plot.

% Write a MATLAB function to make this P-value plot. 

% Remove duplicates and sort p-values from large to small
p = sort(unique(pValues),1,'descend'); 
% Create an array of NaN's to hold Np
Np = ones(length(p),1) .* NaN;

% A loop over all unique p-values to calulate the sum of the p-values that are
% larger than the kth pvalue
for k = 1: length(p)
    Np(k) = sum(pValues > p(k));
end
figure, plot(1-p,Np,'k.'); hold on;
xlabel('1-p'); ylabel('N_p');
title('P-value plot according to Schweder & Spjotvoll 1982');

%% Step 5. P-value plots of REAL microarray data

load prostatecancerexpdata
[H,pValues] = ttest2(independentData', dependentData'); % No 'for' loop!
[T0,T0sd] = PvalPlot(pValues',1);

%% Step 6. Identify significantly modulated genes and plot and label them.

% The slope of our p-value plot, T0, is our estimate for the true number of
% null hypotheses in our data set. However, as S&S point out, there is some
% uncertainty in this effort, largely due to the correlations between
% p-values. In the 1982 paper, they develop a formula for estimating the
% variance in T0, which is returned as 'T0sd' by PvalPlot. One might also
% consider using simulations to estimate the variance of T0 under H0 (Good
% future exercise!)

% Let's take as our estimate of the true # of H0s as T0 + 2*T0sd, meaning
% that the true number of significant differences is:
nSigTrue = size(dependentData,1) - round(T0 + 2*T0sd);

% Now we want to identify the corresponding probes:
[sortedPvals,PvalIndices] = sort(pValues);
% The default for 'sort' is ascending order, so the significant p-values
% are the first ones. We don't really care about the p-values themselves at
% this point; we want to use the indices to find the names and the data so
% that we can plot them:
figure;
% first plot all of the data in black . . .
plot(mean(independentData,2), mean(dependentData,2), 'k.'); hold on;
% then our significant data in red
plot(mean(independentData(PvalIndices([1:nSigTrue]),:),2), ...
     mean(dependentData(PvalIndices([1:nSigTrue]),:),2),'r.');
 
% Find the biggest difference and label that point:
allDiffs = abs(mean(dependentData,2) - mean(independentData,2));
[maxDiff,iMax] = max(allDiffs);
text(mean(independentData(iMax,:)), mean(dependentData(iMax,:)), ...
    probesetIDs(iMax));

% Why does the 'a' look funny? (Ans.: It's a LaTeX thing.)
% Why is this point not colored red? (Hint: Look at the data.)

%% Step 7. Other tests: false discovery rates and "q"

% The beautiful thing about the P-value plot is that it let's us estimate
% the rate of truly null features, based on the common sense notion that
% values to the left on the P-value plot (i.e. high P-values) consist
% almost exclusively of true H0s and the fact that P is distributed
% uniformly under H0. Once we have this estimate, we can calculate the
% "false discovery rate," or "q". While there are various different methods
% to do this, they all rely on the fundamental insights of the P-value
% plot. If you are interested in pursuing this, below are the references
% behind the FDR tests performed in MATLAB's Bioinformatics Toolbox.

figure;
[fdr, q] = mafdr(pValues, 'showplot', true); %Estimate false discovery rate (FDR) of differentially expressed genes from two experimental conditions or phenotypes

% Storey, J.D. (2002). A direct approach to false discovery rates.
% Journal of the Royal Statistical Society 64(3), 479–498.
% 
% *Storey, J.D., and Tibshirani, R. (2003). Statistical significance for
% genomewide studies. Proc Nat Acad Sci 100(16), 9440–9445.
% 
% Storey, J.D., Taylor, J.E., and Siegmund, D. (2004). Strong control
% conservative point estimation and simultaneous conservative consistency
% of false discovery rates: A unified approach. Journal of the Royal
% Statistical Society 66, 187–205.
% 
% Benjamini, Y., and Hochberg, Y. (1995). Controlling the false
% discovery rate: A practical and powerful approach to multiple testing.
% Journal of the Royal Statistical Society 57, 289–300.