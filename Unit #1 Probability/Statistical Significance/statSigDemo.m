% statSigDemo.m
%
% Noodling around with some ideas for an exercise on statistical
% significance; a common error and some intuitions
%
% RTB wrote it, 01 December 2017, home with a head cold

%% Gelman & Stern 2006, The American Statistician (2006) 60:328-331

% The Difference Between “Significant” and “Not Significant” is not
% Itself Statistically Significant
% 
% "Consider two independent studies with effect estimates and standard
% errors of 25 ± 10 and 10 ± 10. The first study is statistically
% significant at the 1% level, and the second is not at all statistically
% significant, being only one standard error away from 0. Thus, it would be
% tempting to conclude that there is a large difference between the two
% studies. In fact, however, the difference is not even close to being
% statistically significant: the estimated difference is 15, with a
% standard error of ?10^2 + 10^2 = 14."
mu1 = 25; mu2 = 10; sem1 = 10; sem2 = 10;
zVal1 = mu1 / sem1;
zVal2 = mu2 / sem2;

% First, how do we think about significance given these numbers? Remember
% that the standard error of the mean is just the standard deviation of the
% sampling distribution for the mean. Because of the Central Limit Theorem,
% we know that the sampling distribution of the mean is Gaussian, so we can
% use our old friends 'normcdf' and 'norminv' to go back and forth between
% probability and standard deviates (i.e. z values).

% So, for the first study, we can see that our mean is 2.5 standard
% deviations from the null value (0), which we just know is more than the
% 1.96 standard deviates required for significance at alpha = 0.05. But how
% do we compute the actual p-value? Recall that the standard normal
% distribution is symetrical, centered around 0 and with s.d. = 1. So we
% can start by asking how much probability mass of N is below -2.5. But
% since we are interested in both tails, we need to multiply this
% probability by 2.
pVal1 = 2 * (normcdf(-zVal1,0,1));    % p = 0.0124

% Note that we could also have used:
% pVal1 = 2 * (1 - normcdf(zVal1,0,1));

% A more literal way to write this, using a graph:
x = -3.5:0.01:3.5;
y = normpdf(x);
figure, plot(x,y,'k-');
hold on
xlabel('z value'); ylabel('Probability');
ax = axis;
line([-zVal1,-zVal1],[ax(3),ax(4)],'Color','r','LineStyle','--');
line([zVal1,zVal1],[ax(3),ax(4)],'Color','r','LineStyle','--');

pRightTail = (1 - normcdf(zVal1,0,1));
pLeftTail = normcdf(-zVal1,0,1);
pVal1 = pRightTail + pLeftTail;

% Now for our 2nd study: mean = 10, sem = 10
pVal2 = 2 * (normcdf(-zVal2,0,1));    % p = 0.3173
line([-zVal2,-zVal2],[ax(3),ax(4)],'Color','b','LineStyle','--');
line([zVal2,zVal2],[ax(3),ax(4)],'Color','b','LineStyle','--');

% plot the data:
figure, errorbar([mu1,mu2],[sem1,sem2]);
hold on
xlabel('Study #'); ylabel('Mean effect');
ax = axis;
line([ax(1),ax(2)],[0,0],'Color','k','LineStyle','--');
axis([ax(1),ax(2),ax(3)-5,ax(4)+5]);
legend('mean +/- SEM')

% What is the sem of the difference? How do we combine variances from
% independent samples?
semDiff = sqrt(sem1^2 + sem2^2);
muDiff = mu1 - mu2;
pValDiff = 2 * normcdf(-(muDiff/semDiff),0,1);  % p = 0.2888

%% Intuition about p-values and error bars

% This is a question that Peter Park uses. How much overlap between error
% bars corresponds to a significant difference?

% Let's start with our existing values and ask how small the sem of each
% study must be until we get signifcance
sem = 11;
myCrit = 0.01;
pValDiff = 0.5;

while pValDiff > myCrit
    sem = sem - 1;
    CI95 = 1.96*sem;
    
    % plot the data:
    figure, errorbar([mu1,mu2],[CI95,CI95]);
    hold on
    xlabel('Study #'); ylabel('Mean effect');
    legend('mean +/- 95% CI');
    ax = axis;
    axis([ax(1),ax(2),ax(3)-5,ax(4)+5]);
    line([ax(1),ax(2)],[0,0],'Color','k','LineStyle','--');
    line([ax(1),ax(2)],[mu1,mu1],'Color','r','LineStyle','--');
    line([ax(1),ax(2)],[mu2,mu2],'Color','b','LineStyle','--');
        
    % What is the sem of the difference? How do we combine variances from
    % independent samples?
    semDiff = sqrt(sem^2 + sem^2);
    muDiff = mu1 - mu2;
    pValDiff = 2 * normcdf(-(muDiff/semDiff),0,1);
    tStr = sprintf('CI_9_5: %.1f, pVal of difference: %.3f',CI95,pValDiff);
    title(tStr);
end

% Could have the students repeat using SEM error bars