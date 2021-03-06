% compareLSRvsML.m: comparing least squares regression vs. maximum
% likelihood fitting
%
% RTB wrote it, 22 September 2020, Northeast Harbor, ME
% reward for submitting Karger review paper

% # of replicates
nReps = 5;

% ground truth y-intercept and slope:
bTrue = [2,1.5];

x = repmat(1:10,nReps,1);
y = bTrue(1) + bTrue(2) .* x;
e = randn(size(y));

% But let's make the last point have an offset and larger variance:
e(:,end) = (randn(nReps,1) .* 7) + 10;
y = y + e;

% Jitter x a little for better visibility:
jFactor = 0.2;
xJitter = (rand(size(y)) .* jFactor) - jFactor/2;

figure
h1 = plot(x+xJitter,y,'bo');
hold on
xlabel('x');
ylabel('y');

% Fit line with least squares regression. I know what it does, because
% Ruilin and I wrote it:
bLSR = [0,0];
[bLSR(2),bLSR(1)] = l_regression(x(:),y(:),'b',0,0.05);

% Fit line with maximum likelihood. Again, I wrote it:
xMean = mean(x);
yMean = mean(y);
ySEM = std(y) ./ sqrt(nReps);
bML = [0,0];
[bML(2),bML(1)] = linfit([xMean',yMean',ySEM']);

% plot the 2 regression lines:
yPredLSR = bLSR(1) + bLSR(2) .* xMean;
h2 = plot(xMean,yPredLSR,'b-');

yPredML = bML(1) + bML(2) .* xMean;
h3 = plot(xMean,yPredML,'r-');

% plot the means of y:
h4 = plot(xMean,yMean,'bo');
set(h4,'MarkerFaceColor','b');

% Does 'fitglm' use ML or LSR?
mdl1 = fitglm(x(:),y(:),'linear','Distribution','normal','Link','identity');
bGLM(1) = mdl1.Coefficients{1,'Estimate'};
bGLM(2) = mdl1.Coefficients{2,'Estimate'};

legend([h1(1),h4,h2,h3],{'raw','mean','LSR fit','ML fit'},'Location','Northwest');

% Print out the different coefficients in a nice table:
T = table([bTrue(1);bLSR(1);bML(1);bGLM(1)],[bTrue(2);bLSR(2);bML(2);bGLM(2)],...
    'VariableNames',{'Beta0','Beta1'},'RowNames',{'True','LSR','ML','GLM'});
display(T);
