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
figure, histogram(manFirst,xBins);
hold on
histogram(womanFirst,xBins);
xlabel({'Percentage of questions asked by women minus';...
    'percentage of attendees who were women (% points)'});
ylabel('Number of seminars');
legend('Man 1st','Woman 1st');
set(gca,'YGrid','on')
title('University seminars, relative share of questions asked by women')

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
% Under what hypothesis should we perform our simulation? Ans.: H0
% What is H0? Ans.: M/W equally likely to ask question regardless of who
% asks the 1st question.
% How do we simlulate this?

% Eventually, I would probably convert the script to a function, in which
% case, the three variables below would be passed to the function:

% 1. Number of questions asked at each seminar:
nQperSeminar = 3;
% 2. Number of seminars to simulate (249 in original study)
nSeminars = length(manFirst) + length(womanFirst);
% 3. Number of simulations to run. When I am first coding up a simulation,
% I generally set this number to something smallish, like 1000. Then, when
% I am confident that it is running correctly, I increase the number to
% 10,000 or 100,000.
nSims = 1000;  

% Variables to hold the results of our simulations:
allEffectSizes = zeros(nSims,1);

% QUESTION (Q1): What is the effect size the authors actually obtained?
realEffectSize = ;

% Setting the random number generator to 'default' ensures that we will all
% get the exact same answer (provided we run the simulation the same number
% of times). You would ordinarily NOT do this. Why?
rng('default');

% TODO: Simulate the experiment nSims times under H0 
for k = 1:nSims
    % Simulate data for one experiment under H0: each row is a question,
    % each column a seminar. Assume a value of '1' means a woman asked the
    % question; '0' means a man asked it.
    !!! Your code here
    
    % Sort according to who asked the 1st questions.
    !!! Your code here
     
    % Now calculate our effect size metric. Assume that attendance is
    % 50/50. So we want to know the proportion of 1's (woman-asked
    % questions) in each column
    !!! Your code here
    
    % Store the effect size for this simulation
    allEffectSizes(k) = ;
end

% plot the results of our simlulation:
figure, histogram(allEffectSizes);
hold on
ax = axis;

% draw a solid black line for the actual effect size
line([realEffectSize,realEffectSize],[ax(3),ax(4)],'Color','k','LineWidth',2);

xlabel('Effect size due to sorting bias');
ylabel('# of simulations');
tStr = sprintf('# of seminars: %d; Questions per seminar: %d',nSeminars,nQperSeminar);
% NOTE: I am an old 'C' programmer, so I like 'sprintf'. But you could
% generate the appropriate text string in a more MATLAB-y way with:
%tStr = ['# of seminars: ' num2str(nSeminars) '; Questions per seminar: ' num2str(nQperSeminar)];
title(tStr);

% QUESTION (Q2): What is the mean effect size, rounded to the nearest whole
% number, if there are only 3 questions asked at each seminar?

% QUESTION (Q3): What if there are 4 questions per seminar?

% QUESTION (Q4): What if there are 5 questions per seminar?

% QUESTION (Q5): Do you see a trend? Try smaller values for nQperSeminar and
% think about what is going on.

% QUESTION (Q6): So then why do the simulation? What does it add?

% QUESTION (Q7):What would happen to our distribution of simulated effect sizes if there
% were 127 seminars instead of 249? Simulate it and compare!

% QUESTION (Q8): What is the smallest number of questions per seminar for
% which you would be 95% confident that the actual effect size obtained
% (25) was not purely do to a sorting bias? Hint: For a given nQperSeminar,
% look at the distribution of your simulated values, and see where the
% actual effect size obtained falls within this distribution.

% QUESTION (Q9): How could you eliminate this bias?

