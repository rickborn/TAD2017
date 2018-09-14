% seminarQuestionBias.m
%
% Do women ask more questions at a seminar if a woman asks the 1st question?
% 
% Prompted by a tweet from Duncan Green referring to a post on his blog:
% http://oxfamblogs.org/fp2p/how-to-stop-men-asking-all-the-questions-in-seminars-its-really-easy/
%
% Original source for data:
% "Women’s visibility in academic seminars: women ask fewer questions than
% men," Alecia Carter, Alyssa Croft, Dieter Lukas, Gillian Sandstrom
%
% RTB wrote it, 14 December 2017
% RTB revised for TAD answers, original name was sexQuestionDemo.m

% Concepts covered:
% 1. histograms for summarizing data
% 2. difference between 'sig' and 'not sig' is not itself nec. stat. sig.
% 3. d-prime as a measure of effect size
% 4. simulating data under H0
% 5. detecting and removing sorting biases
% 6. extracting p-values from simulations

%% Load and plot data

% Each datum represents a value derived from one academic seminar. Values
% are percentage of questions from women minus the percentage of seminar
% attendees who were women. Positive values indicate that women asked more
% questions; negative values indicate that men asked more questions. The
% two variables are the values for each seminar when a woman asked the 1st
% question ('womanFirst') vs. when a man asked the 1st question
% ('manFirst').

% NOTE: I reverse-engineered the raw data based on the graphic in the piece
% by The Economist, so the numbers might not be exactly correct
%
% womanFirst = [60,44,36,36,28,28,28,28,28,24,20,16,16,12,8,8,8,8,8,8,8,8,...
%     4,4,4,4,4,4,0,0,0,0,0,0,0,-4,-4,-4,-4,-4,-4,-8,-8,-8,-8,-8,...
%     -12,-12,-12,-12,-12,-12,-12,-12,-12,-16,-16,-16,-16,-16,-16,-16,-16,-16,...
%     -20,-20,-20,-20,-20,-20,-20,-20,-24,-24,-28,-32,-36,-40];
% 
% manFirst = [24,24,20,20,16,16,16,16,16,16,16,12,8,4,4,4,0,0,0,0,0,...
%     -4,-4,-4,-4,-4,-4,-4,-4,-4,-8,-8,-8,-8,-8,-8,-12,-12,-12,-12,-12,...
%     -12,-12,-12,-12,-12,-12,-12,-12,-12,-12,-16,-16,-16,-16,-20,-20,...
%     -20,-20,-20,-20,-20,-24,-24,-24,-24,-24,-24,-24,-24,-24,-24,-24,-24,...
%     -24,-24,-24,-24,-24,-24,-24,-24,-24,-24,-28,-28,-28,-28,-28,-28,-28,-28,-28,...
%     -32,-32,-32,-32,-32,-32,-32,-32,-32,-32,-36,-36,-36,-36,-36,-36,-36,...
%     -36,-36,-36,-40,-40,-40,-40,-40,-40,-40,-40,-40,-40,-40,-44,-44,-44,-44,...
%     -44,-44,-44,-44,-44,-44,-44,-44,-48,-48,-48,-48,-48,-48,-48,-48,-52,-52,...
%     -52,-52,-52,-52,-52,-56,-56,-56,-56,-56,-56,-56,-56,-60,-60,-60,-60,...
%     -64,-64,-64,-64,-68,-72,-72,-76];

load seminarQuestionData
allSeminars = [womanFirst,manFirst];

% a few useful numbers:
nWF = length(womanFirst);
nMF = length(manFirst);
nTotal = nWF + nMF;

xBins = min(allSeminars):4:max(allSeminars);
figure
subplot(3,1,1);
histogram(manFirst,xBins);
hold on
histogram(womanFirst,xBins);
xlabel({'Percentage of questions asked by women minus';...
    'percentage of attendees who were women (% points)'});
ylabel('Number of seminars');
legend('Man 1st','Woman 1st');
set(gca,'YGrid','on')
title('University seminars, relative share of questions asked by women')

%% Calculate an effect size, d-prime

% NOTE: You were not asked to do this in class, but it is a more common
% measure of effect size than the simple difference in the means (which is
% what we used. It is just the raw difference divided by a pooled estimate
% of the standard deviation.

% NOTE: This formula assumes that the two distributions have the same
% variance. We can test for this with an F-test:
[h,p] = vartest2(womanFirst,manFirst);  % p = 0.1123
% We fail to reject H0 of equal variances, so we can proceed as if they are equal.

% https://en.wikipedia.org/wiki/Sensitivity_index
dPrimeReal = (mean(womanFirst) - mean(manFirst)) / ...
    sqrt(0.5*(var(womanFirst) + var(manFirst)));

tStr = sprintf('dPrime = %0.2f',dPrimeReal);
ax = axis;
text(-70,2/3*ax(4),tStr);

%% Do women tend to ask more questions when a woman asks the first question?

% NOTE: Here is the original tweet from Duncan Green:
% In academic seminars, ‘Men are > 2.5 times more likely to pose questions
% to the speakers. This male skew was observable only in those seminars in
% which a man asked first question. When a woman did so, gender split
% disappeared’. CHAIRS PLEASE NOTE – FIRST Q TO A WOMAN – EVERY TIME.’

% It was based on this story he had read in The Economist:
%
% (https://www.economist.com/news/science-and-technology/ ...
% 21732082-there-easy-fix-women-ask-fewer-questions-men-seminars)
% Print edition | Science and technology
% December 7, 2017
%
% "ONE theory to explain the low share of women in senior academic jobs is
% that they have less self-confidence than men. This hypothesis is
% supported by data in a new working paper, by a team of researchers from
% five universities in America and Europe. In this study, observers counted
% the attendees, and the questions they asked, at 247 departmental talks
% and seminars in biology, psychology and philosophy that took place at 35
% universities in ten countries. On average, half of each seminar’s
% audience was female. Men, however, were over 2.5 times more likely to
% pose questions to the speakers—an action that may be viewed (rightly or
% wrongly) as a sign of greater competence.
% 
% This male skew in question-asking was observable, however, only in those
% seminars in which a man asked the first question. When a woman did so,
% the gender split in question-asking was, on average, proportional to that
% of the audience. Simply handing the microphone to a woman rather than a
% man when the floor is opened for questions may make a difference, however
% small, to one of academia’s most intractable problems."

% This is a classic example of the fallacy that Gelman & Stern write about:
% The Difference Between “Significant” and “Not Significant” is not Itself
% Statistically Significant
%
% In other words, what the blogger did was two t-tests, one of which showed
% that the man-first data are highly statistically significant from 0. That
% is, men ask significantly more questions when a man asks the 1st Q.
[~,pMF] = ttest(manFirst);  % p = 3.56e-36

% But this effect "disappeared" (i.e. was not sig. diff. from 0) when a
% woman asked the 1st question:
[~,pWF] = ttest(womanFirst);  % p = 0.31

% But this does not necessarily mean that the two groups are significantly
% different FROM EACH OTHER. This is very common statistical error that is
% beautifully described in a classic paper:
%
% Gelman A & Stern H (2006) "The Difference Between 'Significant' and 'Not
% Significant' is not Itself Statistically Significant", The American
% Statistician (2006) 60:328-331


% To establish a difference between the two groups, you need to directly
% compare them. A simple way is with a 2-sample t-test
[~,pDiff,~,stats] = ttest2(womanFirst,manFirst);    % p = 4.6e-16

% The above test assumes equal variances. We can test for this directly
% using a two-sample F-test for equal variances:
[hVar,pVar] = vartest2(womanFirst,manFirst);

% Note that the t-statistic (returned in 'stats.tstat') is also a kind of
% effect size. It is normalized by a pooled estimate of the s.e.m., so it
% will be larger than our d-prime, which is normalized by the pooled
% estimate of the s.d.

% or, if we suspect the data are not normally distributed, we can use a
% non-parametric test based on the ranks of the data. This is known as the
% Wilcoxon Rank-Sum Test (equivalent to the Mann-Whitney U Test):
pRankSum = ranksum(manFirst,womanFirst);    % p = 5.9e-15

% In any case, there does appear to be a very real effect of who asks the
% first question. So the conclusion is correct, even though it was
% originally based on a faulty comparison.

%% Bias produced by the way the data were sorted?

% Finally, note that there is a possible source of bias in the way the data
% have been stratified (a fancy statistical word for "sorted"). If there
% are a relatively small number of questions asked at any given seminar,
% then when you divide the data into two subsets where the sex of the 1st
% question asker is fixed, you create a bias in that direction. We can
% simlulate this to see how big the bias effect is.

% Keys to the simulation:
% Under what hypothesis should we peroform our simulation? Ans.: H0
% What is H0? Ans.: M/W equally likely to ask question regardless of who
% asks the 1st question.
% How do we simlulate this?

% Eventually, I would probably convert the script to a function, in which
% case, the below would be variables passed to the function:
nQperSeminar = 6;
nSeminars = length(manFirst) + length(womanFirst); % i.e. 249 in original study
nSims = 10000;

% Variables to hold the results of our simulations:
allEffectSizes = zeros(nSims,1);
allDPrimes = zeros(nSims,1);

% Calculate the effect size the authors actually obtained:
realEffectSize = mean(womanFirst) - mean(manFirst);

% Setting the random number generator to 'default' ensures that we will all
% get the exact same answer (provided we run the simulation the same number
% of times). You would ordinarily NOT do this. Why?
rng('default');
for k = 1:nSims
    % Simulate data: each row is a question, each column a seminar. Assume a
    % value of '1' means a woman asked the question; '0' means a man asked it.
    allData = round(rand(nQperSeminar,nSeminars));
    
    % Sort according to who asked the 1st questions. We can simulate either
    % with the original error (i.e. counting all rows) or with the 'fix',
    % which is just to exclude the 1st question.
    fixFlag = 0;    % fix the bias (i.e. make it go away)
    if fixFlag
        simManFirst = allData(2:end,allData(1,:) == 0);
        simWomanFirst = allData(2:end,allData(1,:) == 1);
    else
        simManFirst = allData(:,allData(1,:) == 0);
        simWomanFirst = allData(:,allData(1,:) == 1);
    end
    
    % Now calculate our metric. For now, assume attendance is 50/50. So we want
    % to know the proportion of 1's (woman-asked questions) in each column
    % NOTE: In the 'sum' commands below, we specify the 1st dimension even
    % though this is the default, because otherwise it will fail when
    % nQperSeminar is 1 and we have a row vector. That is, for an m-by-1
    % vector, a=[1,2,3,4], 'sum(a)' will return 10, whereas 'sum(a,1)' will
    % return [1 2 3 4]. This is a teachable moment.
    simMFscores = round(((sum(simManFirst,1) ./ nQperSeminar) - 0.5) * 100);
    simWFscores = round(((sum(simWomanFirst,1) ./ nQperSeminar) - 0.5) * 100);
    
    allEffectSizes(k) = mean(simWFscores) - mean(simMFscores);
    allDPrimes(k) = (mean(simWFscores) - mean(simMFscores)) / ...
        sqrt(0.5*(var(simWFscores) + var(simMFscores)));
end

% plot the results of our simlulation:
subplot(3,1,2);
histogram(allEffectSizes);
hold on
ax = axis;

% draw a solid black line for the actual effect size
line([realEffectSize,realEffectSize],[ax(3),ax(4)],'Color','k','LineWidth',2);

% Calculate a 95% confidence interval for the simulation:
myAlpha = 0.05;     % MATLAB convention for determining 95% CI
idxHi = ceil(nSims * (1 - myAlpha/2));
idxLo = floor(nSims * (myAlpha/2));
sortedEffectSizes = sort(allEffectSizes);
simCI = [sortedEffectSizes(idxLo),sortedEffectSizes(idxHi)];
line([simCI(1),simCI(1)],[ax(3),ax(4)],'Color','b','LineStyle','--');
line([simCI(2),simCI(2)],[ax(3),ax(4)],'Color','b','LineStyle','--');

xlabel('Effect size');
ylabel('# of simulations');
tStr = sprintf('# of seminars: %d; Questions per seminar: %d',nSeminars,nQperSeminar);
% NOTE: I am an old 'C' programmer, so I like 'sprintf'. But you could
% generate the appropriate text string in a more MATLAB-y way with:
%tStr = ['# of seminars: ' num2str(nSeminars) '; Questions per seminar: ' num2str(nQperSeminar)];
title(tStr);

pValue1 = sum(allEffectSizes >= realEffectSize) / nSims;
% Why do I do this? Can a p-value ever be 0?
if pValue1 == 0
    pValue1 = 1 / (nSims + 1);
end

% print p-value on plot
xTxt = (floor(ax(1)/5) + 1) * 5;
yTxt = 0.75 * ax(4);
txtStr = sprintf('p = %0.3f',pValue1);
text(xTxt,yTxt,txtStr);

% QUESTION: What is the mean effect size, rounded to the nearest whole
% number, if there are only 4 questions asked at each seminar?
% ANSWER: 25

% QUESTION: What if there are 5 questions per seminar?
% ANSWER: 20

% QUESTION: What if there are 6 questions per seminar?
% ANSWER: 17

% QUESTION: Do you see a trend? Try smaller values for nQperSeminar and
% think about what is going on.
% ANSWER: We can reason it out that, if nQperSeminar is 1, then all values
% must be 100. For nQperSeminar = 2, on average we should get a value of
% 50, and so on: meanBiasEffect = 100 / nQuestionsPerSeminar

% QUESTION: So then why do the simulation?
% ANSWER: By simulating the conditions we used in our actual experiment
% (i.e. nSeminars = 249), we get an idea of the random variability in the
% expected effect size under H0.

% What would happen to our distribution of simulated effect sizes if there
% were 127 seminars instead of 249?
% ANSWER: The variance would increase.

% QUESTION: What is the smallest number of questions per seminar for which
% you would be 95% confident that the actual effect size obtained (25) was
% not purely do to a sorting bias?
% ANSWER: 5

% QUESTION: How could you eliminate this bias?
% ANSWER1: Sort by the 1st question, but exclude it from the calculation.
% ANSWER2: Another approach might be to subtract the mean bias effect. But
% this would not take into account the variability due to sampling error.
%% Bonus: Histogram of simulated d-primes

subplot(3,1,3);
histogram(allDPrimes);
hold on
ax = axis;
line([dPrimeReal,dPrimeReal],[ax(3),ax(4)],'Color','k','LineWidth',2);
xlabel('d-prime');
ylabel('# of simulations');
tStr = sprintf('Questions per seminar = %d', nQperSeminar);
title(tStr);

pValue2 = sum(allDPrimes >= dPrimeReal) / nSims;
if pValue2 == 0
    pValue2 = 1 / (nSims + 1);
end

% print p-value on plot
xTxt = (0.1 * (ax(2) - ax(1))) + ax(1);
yTxt = 0.75 * ax(4);
txtStr = sprintf('p = %0.2f',pValue2);
text(xTxt,yTxt,txtStr);

%% Bonus: Variable number of questions per seminar

% Let's suppose we don't want to assume that there were the exact same
% number of questions asked at each seminar. How could we build this into
% our simulation?

% Maybe we just know that, on average, there were m questions per seminar.
% Since we can assume that seminars have no memory (i.e. for how many
% questions were asked at the previous seminars), we can just draw values
% from a Poisson distribution with lambda set to the mean number of
% questions over all seminars. For each round of our simulation, we will
% first use 'Poissrnd' to generate the number of questions asked at each of
% the 249 seminars. Then we will generate a matrix big enough to
% accommodate the largest number of questions, then replace data values
% with NaN's to tailor each seminar to the actual number of questions. In
% MATLAB, 'NaN' stands for 'Not a Number' and is useful for representing
% missing data. Most arithmetical functions have graceful ways of dealing
% with NaN's such that they are not counted. See code below for how this is
% done.

% Variables to hold the results of our simulations:
allEffectSizes = zeros(nSims,1);
allDPrimes = zeros(nSims,1);

% Calculate the effect size the authors actually obtained:
realEffectSize = mean(womanFirst) - mean(manFirst);

rng('default');
for k = 1:nSims
    % How many questions were asked at each seminar?
    allQperSeminar = poissrnd(nQperSeminar,1,nSeminars);
    maxQ = max(allQperSeminar);
    
    % Simulate data: each row is a question, each column a seminar. Assume a
    % value of '1' means a woman asked the question; '0' means a man asked it.
    allData = round(rand(maxQ,nSeminars));
    
    % Now trim each seminar to the correct number using NaN's. It doesn't
    % matter where we trim, as long as we don't mess with the first row, so
    % we'll put them at the end:
    for j = 1:nSeminars
        % How many questions do we need to replace in the seminar (column)?
        nReplace = maxQ - allQperSeminar(1,j);
        % Create a column of NaN's of the appropriate size:
        thisNaNpad = ones(nReplace,1) .* NaN;
        % Paste this over the end of the data column, sparing the appropriate
        % number of questions (rows) above:
        allData(allQperSeminar(1,j)+1:end,j) = thisNaNpad;
    end
    
    % Sort according to who asked the 1st questions. We can simulate either
    % with the original error (i.e. counting all rows) or with the 'fix',
    % which is just to exclude the 1st question.
    fixFlag = 0;    % fix the bias (i.e. make it go away)
    if fixFlag
        simManFirst = allData(2:end,allData(1,:) == 0);
        simWomanFirst = allData(2:end,allData(1,:) == 1);
    else
        simManFirst = allData(:,allData(1,:) == 0);
        simWomanFirst = allData(:,allData(1,:) == 1);
    end
    
    % Now calculate our metric. For now, assume attendance is 50/50. So we want
    % to know the proportion of 1's (woman-asked questions) in each column
    % NOTE: We need to make two changes to our previous code. First, we
    % include 'omitnan' as an argument passed to the 'sum' command so that
    % it will ignore the 'NaN' values. Second, we need to divide each
    % column sum by the actual number of questions simulated in that colum.
    % To do that we just count up the number of values in each column that
    % are not NaN's.
    simMFscores = round(((sum(simManFirst,1,'omitnan') ./ sum(~isnan(simManFirst))) - 0.5) * 100);
    simWFscores = round(((sum(simWomanFirst,1,'omitnan') ./ sum(~isnan(simWomanFirst))) - 0.5) * 100);
    
    allEffectSizes(k) = mean(simWFscores) - mean(simMFscores);
    allDPrimes(k) = (mean(simWFscores) - mean(simMFscores)) / ...
        sqrt(0.5*(var(simWFscores) + var(simMFscores)));
end

% plot the results of our simlulation:
figure
subplot(3,1,1);
histogram(allEffectSizes);
hold on
ax = axis;

% draw a solid black line for the actual effect size
line([realEffectSize,realEffectSize],[ax(3),ax(4)],'Color','k','LineWidth',2);

% Calculate a 95% confidence interval for the simulation:
myAlpha = 0.05;     % MATLAB convention for determining 95% CI
idxHi = ceil(nSims * (1 - myAlpha/2));
idxLo = floor(nSims * (myAlpha/2));
sortedEffectSizes = sort(allEffectSizes);
simCI = [sortedEffectSizes(idxLo),sortedEffectSizes(idxHi)];
line([simCI(1),simCI(1)],[ax(3),ax(4)],'Color','b','LineStyle','--');
line([simCI(2),simCI(2)],[ax(3),ax(4)],'Color','b','LineStyle','--');

xlabel('Effect size');
ylabel('# of simulations');
tStr = sprintf('# of seminars: %d; Avg. # Q per seminar: %d',nSeminars,nQperSeminar);
title(tStr);

pValue3 = sum(allEffectSizes >= realEffectSize) / nSims;
% Why do I do this? Can a p-value ever be 0?
if pValue3 == 0
    pValue3 = 1 / (nSims + 1);
end

% print p-value on plot
xTxt = (floor(ax(1)/5) + 1) * 5;
yTxt = 0.75 * ax(4);
txtStr = sprintf('p = %0.3f',pValue3);
text(xTxt,yTxt,txtStr);

% This seems to have a very profound effect on our simulation. Did we make
% a mistake? To start to get some insight into this question, let's look at
% the distribution of the actual number of questions asked:
subplot(3,1,2);
histogram(allQperSeminar);
xlabel('# of Questions asked');
ylabel('# of Seminars');

% This is the Poisson distribution for a lambda of 5. Right away we can see
% that there is more mass for the smaller numbers of questions:
% sum(allQperSeminar <= 5) is about 147 (will vary per simulation)
% sum(allQperSeminar > 5) is about 102
% So this pushes our biases towards bigger values than when we assumed that
% all seminars had exactly 5 seminars.

% But there is something else going on, too. Recall the formula we derived
% for the mean bias as a function of the number of questions asked per
% seminar:
%
% meanBiasEffect = 100 / nQuestionsPerSeminar
%
% Because the # of questions asked is in the denominator, it means that
% seminars in which very few questions were asked produce a HUGE bias,
% while those with a large number of questions produce a progressively
% smaller bias. This also skews the overall distribution towards larger
% biases.

% Plot the effect of # Q per seminar on the mean bias:
xQuestions = 1:maxQ;
yBias = 100 ./ xQuestions;
subplot(3,1,3);
plot(xQuestions,yBias,'b-','LineWidth',2);
xlabel('# of Questions asked');
ylabel('Mean bias');

% I would never have come to these conclusions without having done the
% simulations! In fact, my intution going in was that the mean of the
% simulated bias distribution would be the same whether we assumed that ALL
% seminars had exactly the same number of questions or whether we assumed
% that they had the same average number that varied randomly from seminar
% to seminar.