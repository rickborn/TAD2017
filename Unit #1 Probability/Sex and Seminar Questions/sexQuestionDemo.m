% sexQuestionDemo.m
%
% Do women ask more questions at a seminar if a woman asks the 1st
% question?
% 
% Prompted by a tweet from Duncan Green referring to a post on his blog:
% http://oxfamblogs.org/fp2p/how-to-stop-men-asking-all-the-questions-in-seminars-its-really-easy/
%
% Original source for data:
% "Women’s visibility in academic seminars: women ask fewer questions than
% men," Alecia Carter, Alyssa Croft, Dieter Lukas, Gillian Sandstrom
%
% RTB wrote it, 14 December 2017

% Each datum represents a value derived from one academic seminar. Values
% are percentage of questions from women minus the percentage of seminar
% attendees who were women. Positive values indicate that women asked more
% questions; negative values indicate that men asked more questions. The
% two variables are the values for each seminar when a woman asked the 1st
% question ('womanFirst') vs. when a man asked the 1st question
% ('manFirst').

%% Load and plot data

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
% different from each other. To establish this, you need to directly
% compare them. A simple way is with a 2-sample t-test
[~,pDiff,~,stats] = ttest2(womanFirst,manFirst);    % p = 4.6e-16

% Note that the t-statistic (returned in 'stats.tstat') is also a kind of
% effect size. It is normalized by a pooled estimate of the s.e.m., so it
% will be larger than our d-prime, which is normalized by the pooled
% estimate of the s.d.

% or, if we suspect the data are not normally distributed
pRankSum = ranksum(manFirst,womanFirst);    % p = 5.9e-15

%% Or we can do a permutation test
realMedDiff = median(womanFirst) - median(manFirst);

nPerm = 100000;
permMedDiff = zeros(nPerm,1);

for k = 1:nPerm
    shuffledData = allSeminars(randperm(nTotal));
    permMedDiff(k) = median(shuffledData(1:nWF)) - median(shuffledData(nWF+1:end));
end

figure
histogram(permMedDiff);
xlabel('Difference in permuted medians');
ylabel('Number');
ax = axis;
line([realMedDiff,realMedDiff],[ax(3),ax(4)],'Color','k','LineStyle','--');

pVal = sum(permMedDiff >= realMedDiff) / nPerm;
if pVal == 0
    pVal = 1/(nPerm+1);
end

% In any case, there does appear to be a very real effect of who asks the
% first question. So the conclusion is correct, even though it was
% originally based on a faulty comparison.

%% Bias produced by the way the data were sorted?

% Finally, note that there is another possible source of bias in the way
% the data have been stratified. If there are a relatively small number of
% questions asked at any given seminar, then when you divide the data into
% two subsets where the sex of the 1st question asker is fixed, you create
% a bias in that direction. We can simlulate this to see how big the bias
% effect is.

% Keys to the simulation:
% Under what hypothesis should we peroform our simulation? Ans.: H0
% What is H0? Ans.: M/W equally likely to ask question regardless of who
% asks the 1st question.
% How do we simlulate this?

nQperSeminar = 5;
nSeminars = length(manFirst) + length(womanFirst); % i.e. 249 in original study
% nSeminars = 127;
nSims = 10000;
allEffectSizes = zeros(nSims,1);
allDPrimes = zeros(nSims,1);

realEffectSize = mean(womanFirst) - mean(manFirst);
% Use ttest2 to get confidence intervals on our real value
[~,~,ci,statsl] = ttest2(womanFirst,manFirst);

rng('default');
for k = 1:nSims
    % Simulate data: each row is a question, each column a seminar. Assume a
    % value of '1' means a man asked the question; '0' means a woman asked it.
    allData = round(rand(nQperSeminar,nSeminars));
    
    % Sort according to who asked the 1st questions. We can simulate either
    % with the original error (i.e. counting all rows) or with the 'fix',
    % which is just to exclude the 1st question.
    fixFlag = 0;    % fix the bias
    if fixFlag
        simManFirst = allData(2:end,allData(1,:) == 1);
        simWomanFirst = allData(2:end,allData(1,:) == 0);
    else
        simManFirst = allData(:,allData(1,:) == 1);
        simWomanFirst = allData(:,allData(1,:) == 0);
    end
    
    % Now calculate our metric. For now, assume attendance is 50/50. So we want
    % to know the proportion of 0's (woman-asked questions) in each column
    % NOTE: In the 'sum' commands below, we specify the 1st dimension even
    % though this is the default, because otherwise it will fail when
    % nQperSeminar is 1 and we have a row vector. That is, for an mx1
    % vector, a=[1,2,3,4], 'sum(a)' will return 10, whereas 'sum(a,1)' will
    % return [1 2 3 4]. This is a teachable moment.
    simMFscores = round(((sum(~simManFirst,1) ./ nQperSeminar) - 0.5) * 100);
    simWFscores = round(((sum(~simWomanFirst,1) ./ nQperSeminar) - 0.5) * 100);
    
    allEffectSizes(k) = mean(simWFscores) - mean(simMFscores);
    allDPrimes(k) = (mean(simWFscores) - mean(simMFscores)) / ...
        sqrt(0.5*(var(simWFscores) + var(simMFscores)));
end

figure
histogram(allEffectSizes);
hold on
ax = axis;

% a solid black line for the actual effect size
line([realEffectSize,realEffectSize],[ax(3),ax(4)],'Color','k','LineWidth',2);

% CI for the real effect size.
% NOTE: This now seems less interesting than the CI for the simulation.
% line([ci(1),ci(1)],[ax(3),ax(4)],'Color','k','LineStyle','--');
% line([ci(2),ci(2)],[ax(3),ax(4)],'Color','k','LineStyle','--');

% CI for the simulation:
myAlpha = 0.05;     % MATLAB convention for determining 95% CI
idxHi = ceil(nSims * (1 - myAlpha/2));
idxLo = floor(nSims * (myAlpha/2));
sortedEffectSizes = sort(allEffectSizes);
simCI = [sortedEffectSizes(idxLo),sortedEffectSizes(idxHi)];
line([simCI(1),simCI(1)],[ax(3),ax(4)],'Color','b','LineStyle','--');
line([simCI(2),simCI(2)],[ax(3),ax(4)],'Color','b','LineStyle','--');

xlabel('Effect size');
ylabel('# of simulations');
tStr = sprintf('# of seminars:%d; Questions per seminar:%d',nSeminars,nQperSeminar);
title(tStr);

pValue1 = sum(allEffectSizes >= realEffectSize) / nSims;
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
% ANSWER2: Another approach might be to subtract the mean bias effect.
%% Histogram of simulated d-primes

figure
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

%% Plot raw simulation data

% figure
% allSimData = [simMFscores,simWFscores];
% 
% xBins = min(allSimData):max(allSimData);
% histogram(simMFscores,xBins);
% hold on
% histogram(simWFscores,xBins);
% xlabel({'Percentage of questions from women minus';...
%     'percentage of attendees who are women (% points)'});
% ylabel('Number of seminars');
% legend('Man 1st','Woman 1st');
% 
% [~,pDiff] = ttest2(simWFscores,simMFscores,'tail','right'); % p = 
% 
% % or, if we suspect the data are not normally distributed
% pRankSum = ranksum(simWFscores,simMFscores); 