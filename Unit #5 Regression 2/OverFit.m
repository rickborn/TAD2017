function [p1,p2,sse1,sse2] = OverFit(nf,po,x)

% OVERFIT.M: Demonstrates problem of over fitting data
%
% Inputs:
% - nf, noise factor (gain for multiplying noise)
% - po, polynomial order to use for over fit
% - x, the independent variable
%
% Outputs:
% - p1, coefficients of regression line
% - p2, coefficients of higher order polynomial
% - sse1, sum of squared errors for poly fit to half of data
% - sse2, sum of squared errors for same poly to other half of data
%
% e.g. [p1,p2,sse1,sse2] = OverFit(10,4,[1:12]);
%
% RTB wrote it, Christmas morning 2010

% Once you've run the demo, see how it extrapolates:
% e.g. compare:
%   polyval(p1, max(x) .* 2)
%   polyval(p2, max(x) .* 2)

if nargin < 3, x = 1:12; end
if nargin < 2, po = 4; end
if nargin < 1, nf = 20; end

% generate linear data with noise and fit with a line
noise = (rand(1,length(x)) .* nf) - (nf/2);
y = x+noise;
plot(x,y,'bo','MarkerFaceColor','b'); hold on;
%[h,p1] = fitline(x,y,1);
[p1,S] = polyfit(x,y,1);
f = polyval(p1,x);
plot(x,f,'r-','LineWidth',2);

% now over fit the same data with higher order polynomial
[p2,S] = polyfit(x,y,po);
f = polyval(p2,x);
plot(x,f,'k-','LineWidth',2);
%keyboard;

% Now see how poorly the overfit model extrapolates:
% e.g. compare:
extrapFactor = 1.5;
y1 = polyval(p1, max(x) .* extrapFactor);    % linear fit
y2 = polyval(p2, max(x) .* extrapFactor);    % over fit
figure, plot(max(x) .* extrapFactor, y1, 'r*');
hold on
plot(max(x) .* extrapFactor, y2, 'k*');
%keyboard;
axis([0 max(x) 0 max(y)]);

% Another way of dealing with this issue is to fit the model to a random
% subset of the data, then see how it fits the other half.
o = randperm(length(x));            % scramble the order of the data
t1 = sort(o(1:floor(length(x)/2)));       % 1st half
t2 = sort(o(floor(length(x)/2)+1:end));   % 2nd half

% fit 1st half
[p3,S] = polyfit(x(t1),y(t1),po);
[f1,e1] = polyval(p3,x(t1),S);
[f2,e2] = polyval(p3,x(t2),S);
[f3] = polyval(p3,x);
figure;
plot(x(t1),y(t1),'bo','MarkerFaceColor','b'); hold on;
plot(x(t1),f1,'b-','LineWidth',2);
%keyboard;
plot(x(t2),y(t2),'ro','MarkerFaceColor','r');
plot(x(t2),f2,'r-','LineWidth',2);
xx = min(x):0.1:max(x);
yy = polyval(p3,xx,S);
plot(xx,yy,'k-','LineWidth',2);
% plot(x,f3,'k-','LineWidth',2);

% calculate sum of squared errors
sse1 = 0; sse2 = 0;
for k = 1:length(t1)
    txx = find(xx == x(t1(k)));
    tx = find(x == x(t1(k)));
    sse1 = sse1 + ((y(tx) - yy(txx)).^2);
end

for k = 1:length(t2)
    txx = find(xx == x(t2(k)));
    tx = find(x == x(t2(k)));
    sse2 = sse2 + ((y(tx) - yy(txx)).^2);
end

% These are not really the residuals--see help for 'polyval' 
% sum(e1)    % sum of residuals for fit data
% sum(e2)    % sum of residuals for the other data