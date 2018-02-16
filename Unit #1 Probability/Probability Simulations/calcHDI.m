function [hdi] = calcHDI(allTH,xBins,myAlpha)

% Find the 95% highest density interval for a sampling distribution:
myDist = hist(allTH,xBins);

% Find the index corresponding to the maximum value in the distribution:
maxPNdx = find(myDist == max(myDist));
% Find the critical value for the area under the curve
critVal = (1 - myAlpha) * trapz(myDist);

% Now we just start in the middle and work our way out symmetrically until
% we surpass 95% of the area:
myVal = 0;
ctr = 0;
while myVal < critVal
    ctr = ctr+1;
    myVal = trapz(myDist(maxPNdx-ctr:maxPNdx+ctr));
end
% HDI = "highest density interval"
hdi = [xBins(maxPNdx - ctr),xBins(maxPNdx + ctr)];


% NOTE: This method is not perfect, since it moves out symmetrically in the
% index values, NOT in area. Insofar as the posterior distributions are not
% symmetrical, the HDIs will be shifted a bit.