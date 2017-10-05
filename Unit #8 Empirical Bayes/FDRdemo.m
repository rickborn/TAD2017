% FDRdemo.m: explores the ideas behind false discovery rates
%
% Based on Biostatistics-Lecture_notes-RTB.docx and ideas in seminal papers
% by Benjamini & Hochberg (1995), Storey (2002) and Storey & Tibshirani
% (2003)
%
% RTB wrote it

% Concepts covere:
% 1. Distribution of p-values under H0
% 2. P-value plot after Schweder & Spjotvoll (1982); 'PvalPlot.m' (RTB)
% 3. Intuition for the numerator and denominator of the FDR calculation
% 4. Benjamini and Hochberg procedure to control FDR
% 5. MATLAB's built-in FDR calculation, 'mafdr', using Storey's approach

% Data example is from Efron (EB) p. 15: Genetic expression levels from N =
% 6033 genes were obtained for n = 102 men, n1 = 50 normal control subjects
% and n2 = 52 prostate cancer patients. The goal of the study was to
% discover a small number of "interesting" genes whose expression levels
% differed between the patients and controls.

%% Load in some real data (courtesy of Brad Efron)
% cd 'C:\usr\rick\doc\Committees\PIN\PIN Director\Courses\Stats\Efron Data'
load prostateData;

%% What does the distribution of p-values look like under H0?
nSims = 100000;
nSamps = 50;
% generate some H0 data
[~,nullPvals] = ttest2(randn(nSamps,nSims),randn(nSamps,nSims));
figure, hist(nullPvals,100);
xlabel('p-value'); ylabel('#');
title('P-values under H0');

% Does the sample size matter?
% Does it matter which test we use?
%% False Discovery Rate (Benjamini & Hochberg 1995): 
%     A	sensitivity: P(T+/D+) vs. positive predictive value (PPV): P(D+/T+)
%     B	FDR is analogous to PPV: P(H0 is true / H0 was rejected)
%     C	Two handy features:
%         1)	FDR is calculated based on the p-values from the multiple tests
%         2)	more appropriate in cases where we expect H0 to be false for a fair number of features (e.g. genetic studies)
%         3)	If we choose an FDR of 0.05, it means that we would like no more than 5% of the features that we classify as significant to be false positives
%               a)	This is not the same as choosing the FWER to be 5%
%               b)	As with FWER, we can either use it to set a criterion or we can use it to adjust our p-values
%     D	FDR_i = (# of false positives) / (# called significant)
%         1)	sort p-values and find the ith ordered p-value (= p_i)
%         2)	if we chose this p-value as the cut-off for significance, we would have i significant p-values
%         3)	if all H0’s are true (uniform distribution of p-values), then the expected # of false positives would be: n x p_i 
%         4)	Thus FDRhat_i = n x p_i / i
%         5)	If all H0 are true, then FDR = FWER
%         6)	FDR is appropriate when a single falsely rejected H0 (i.e. false positive) is not too dangerous

% PvalPlot(nullPvals,1);

% caculate p-values for the real data and plot them
[~,dataPvals] = ttest2(patients',controls');
figure, hist(dataPvals,100);
xlabel('p-value'); ylabel('count');
title('P-values for Prostrate Data');

PvalPlot(dataPvals,1);

%figure, plot(dataPvalsSorted);
% hold on;
% xlabel('rank'); ylabel('p-value');

%% FDR = (# of p-values < crit expected under H0) / (# of p-values < crit in data)
% Given that distribution of p-values under H0 is uniform, the expected #
% of false positives (numerator) is just n x p_i. Work this out:
% Let's just see what happens as we change our criterion
% crit = [0.0005,0.001,0.005, 0.01, 0.05, 0.10];
crit = logspace(-4,-1,10);
H0pValsLessThanCrit = zeros(size(crit));
for k = 1:length(crit)
    H0pValsLessThanCrit(k) = sum(nullPvals < crit(k)) / length(nullPvals);
end
figure, p1 = loglog(crit,H0pValsLessThanCrit,'bo');
hold on;
loglog(crit,H0pValsLessThanCrit,'k-');
xlabel('criterion'); ylabel('proportion of p-values < crit');

%% Now calculate the same thing for the actual data
HApValsLessThanCrit = zeros(size(crit));
for k = 1:length(crit)
    HApValsLessThanCrit(k) = sum(dataPvals < crit(k)) / length(dataPvals);
end
p2 = loglog(crit,HApValsLessThanCrit,'ro');
loglog(crit,HApValsLessThanCrit,'k-');

legend([p1,p2],'H0','Data','Location','NorthWest');
%% The FDR is just the ratio of the two measures or:
FDR = H0pValsLessThanCrit ./ HApValsLessThanCrit;
FDRpercent = round(FDR .* 1000) ./ 10;
text(crit',(HApValsLessThanCrit .* 1.5)', num2str(FDRpercent'), 'HorizontalAlignment','Center');

%% Use the Benjamini and Hochberg procedure to control FDR
% desired FDR = q*
qStar = 0.10;
% let m be the total # of hypotheses tested
nTotal = length(dataPvals);
% sort p-values from smallest to largest
[dataPvalsSorted,~] = sort(dataPvals,'ascend');
% find the largest p-value that is less than i/m*q
for k = 1:nTotal
    if dataPvalsSorted(k) > ((k/nTotal)*qStar)
        break
    else
        continue
    end
end
k = k-1

%% Convince yourself that k=59 is right
fdr58 = (dataPvalsSorted(58) * nTotal) / 58
fdr59 = (dataPvalsSorted(59) * nTotal) / 59
fdr60 = (dataPvalsSorted(60) * nTotal) / 60

%% Plot the first 100 ordered p-values (Efron, EB, p. 50)
iPvals = 1:100;
figure, plot(iPvals,dataPvalsSorted(iPvals),'k.');
hold on
xlabel('index i'); ylabel('P-value');

% calculate the rejection region for qStar = 0.10
pValReject = (qStar .* iPvals) ./ nTotal;
h1 = plot(iPvals,pValReject,'r-');

% calculate the rejection region for qStar = 0.05
q05 = 0.05;
pValReject = (q05 .* iPvals) ./ nTotal;
h2 = plot(iPvals,pValReject,'b-');

legend([h1,h2],{'qStar = 0.10','qStar = 0.05'},'Location','Northwest');

%% Repeat the above, but on a semilog plot
figure, semilogy(iPvals,dataPvalsSorted(iPvals),'k.');
hold on
xlabel('index i'); ylabel('P-value');

% calculate the rejection region for qStar = 0.10
pValReject = (qStar .* iPvals) ./ nTotal;
semilogy(iPvals,pValReject,'r-');

% calculate the rejection region for qStar = 0.05
q05 = 0.05;
pValReject = (q05 .* iPvals) ./ nTotal;
semilogy(iPvals,pValReject,'b-');

%% Alternative algorithm for B & H calculation of FDR
% desired FDR = q*
qStar = 0.10;
% let m be the total # of hypotheses tested
nTotal = length(dataPvals);
% sort p-values from smallest to largest
[dataPvalsSorted,dataPvalIndices] = sort(dataPvals,'ascend');
% find the largest p-value that is less than i/m*q
k = 1;
while dataPvalsSorted(k) <= ((k/nTotal)*qStar)
    k = k + 1;
end
k = k-1
%% Now let's use the FDR functions in MATLAB
% pValues = mattest(patients,controls,'permute',true);
figure
[mlFDR,q,Pi0] = mafdr(dataPvals,'showplot',true);

% Pi0 is the estimated a-priori probability that the null hypothesis is true

%% Let's do a simulation of FDR performance
% We control the proportion of null and non-null z-values
% See fig. 4.4, pp. 55-56 of Efron's book on Empirical Bayes

nNull = 2850; nNonNull = 150;
Pi0true = nNull / (nNull + nNonNull);
nSims = 1000;

allFDP = zeros(nSims,1);
% draw null values from std. normal dist
zNull = normrnd(0,1,nNull,1);
zNonNull = normrnd(2.5,1,nNonNull,1);
zAll = [zNull;zNonNull];    % mix them together

