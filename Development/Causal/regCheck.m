% regCheck.m
%
% Reading Judea Pearl's book on causality, "The Book of Why," and I'm
% confused by the way he appears to be using interchangably the regression
% coefficients (i.e. beta's) and the correlation coefficients.
%
% RTB wrote it, 31 May 2019, North Hero, Lake Champlain (Ed Howard's place)

myBeta = [-1,2]';    % slope, intercepts
X = 0:50;
myX = [ones(1,length(X));X]';

% Two data sets with different amounts of noise
y1 = (myX * myBeta) + (randn(length(X),1) .* 5);
y2 = (myX * myBeta) + (randn(length(X),1) .* 20);

figure;
plot(X,y1,'o',X,y2,'s');
hold on

% plot the 2 regression lines:
lsline

% Now do regression
[b1,bint,r,rint,stats1] = regress(y1,myX);
[b2,bint,r,rint,stats2] = regress(y2,myX);

% The regression recovers roughly the same beta's . . .
display([b1,b2]);

% . . . but very different correlation coefficients (R^2)
display([stats1(1),stats2(1)]);