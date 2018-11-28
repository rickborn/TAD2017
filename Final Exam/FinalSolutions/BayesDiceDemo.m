% BayesDiceDemo.m
%
% RTB wrote it, 09 Oct. 2017, rainy Columbus day

% A brown paper bag contains 5 "Platonic solid" dice, each with a different
% number of sides: 4-sided, 6-sided, 8-sided, 12-sided and 20-sided. The
% sides of each die are numbered 1:n, where n is the number of sides, and
% each face is equally likely to come up when the die is rolled. A Bayesian
% statistician randomly selects one of the dice from the bag and rolls it
% behind a screen. She informs you that the roll was a '5'. 

% Q1: Which of the 5 dice did she most likely choose?
nSides = [4,6,8,12,20];
nDice = 5;
prior1 = ones(1,nDice) ./ nDice;
like1 = [0,1/6,1/8,1/12,1/20];
post1 = prior1 .* like1;
post1Nl = post1 ./ sum(post1);

% A1: 6-sided
maxLike = nSides(post1 == max(post1));

% Q2: What is the probability that the 4-sided die was chosen?
p4sided = post1Nl(nSides == 4);
% A2: 0

% Q3: What is the probability that the 8-sided die was chosen?
p8sided = post1Nl(nSides == 8);
% A3: 0.2941

% The chosen die is rolled a 2nd time, and a value of 10 is obtained:
like2 = [0,0,0,1/12,1/20];
post2 = post1Nl .* like2;
post2Nl = post2 ./ sum(post2);

% Q4: What are the odds in favor of the 12-sided die being the chosen one?
odds12sided = post2Nl(nSides==12) / (1 - post2Nl(nSides==12));
% A4 = 2.7778

% Q5: You are now told that there exists a body of literature on the
% tactile preferences of human beings (including Bayesian statisticians)
% showing a marked preference for choosing smooth objects over pointy
% ones. In particular, for dice in brown paper bags, the probability of
% choosing a given die is directly proportional to its number of sides.
% This means, for example, that the 20-sided die is 5/3 more likely to be
% chosen than the 12-sided die. If you were to now see the same sequence of
% two rolls (5,10) from a non-randomly (i.e. biased by human shape
% preferences) chosen die, what are the new odds in favor of the 12-sided
% die being the chosen one?

% One approach is to calculate a new 'prior1' and just re-do the entire
% calculation:
prior1 = nSides ./ sum(nSides);
% Ans = 1.6667

% But, of course, we could just go back and adjust our original odds by
% multiplying by 3/5:
newOdds12sided = odds12sided * (3/5);    % 1.6667



