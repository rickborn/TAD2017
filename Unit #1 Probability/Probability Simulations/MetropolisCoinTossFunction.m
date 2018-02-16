function [allTH,hdi95] = MetropolisCoinTossFunction(z,N,mySigma,a,b,nSims,pFlag)

% MetropolisCoinTossFunction: samples posterior distribution for obtaining
% z heads on N coin tosses.
%
% [allTH,hdi95] = MetropolisCoinTossFunction(z,N,mySigma,a,b,nSims)
%
% Inputs:
% - z, # of heads obtained on . . .
% - N, # of coin tosses
% - mySigma, sigma of the proposal distribution
% - a, b, parameters of the beta distribution for the prior
% - nSims, # of samples to produce
% - pFlag, 1 = make some plots
%
% Outputs:
% - allTH, complete sampling distribution from the posterior
% - hdi95, the 95% highest density interval
%
% e.g. [allTH,hdi95] = MetropolisCoinTossFunction(14,20,0.2,1,1,10000,1);
%
% RTB wrote it,17 January 2018
% based on Kruschke, pp. 158-161, figure 7.4

rng shuffle

allTH = zeros(nSims,1);
THcur = rand;
allTH(1) = THcur;
nAccepted = 0;      % counter for # of times the proposed jump is accepted

for k = 2:nSims
    % randomly generate a proposed jump by drawing from a normal
    % distribution with mean=0 and sigma=mySigma
    THpro = THcur + normrnd(0,mySigma);
    
    if THpro < 0 || THpro > 1
        pMove = 0;
        
    else
        % calculate non-normaized posterior probabilities
        pTHpro = binopdf(z,N,THpro) * betapdf(THpro,a,b);
        pTHcur = binopdf(z,N,THcur) * betapdf(THcur,a,b);
        pMove = min([1,pTHpro/pTHcur]);
    end
    
    if rand < pMove
        allTH(k) = THpro;
        nAccepted = nAccepted + 1;
    else
        allTH(k) = THcur;
    end
    THcur = allTH(k);
end

% calculate the mode
xBins = 0:0.01:1;
myDist = hist(allTH,xBins);
myMode = xBins(myDist == max(myDist));

% calculate the 95% HDI
myAlpha = 0.05;
allTHsorted = sort(allTH);
idxHi = ceil(nSims * (1 - myAlpha/2));
idxLo = floor(nSims * (myAlpha/2));
hdi95 = [allTHsorted(idxLo),allTHsorted(idxHi)];

if pFlag
    figure, subplot(3,1,1);
    xBins = 0:0.01:1;
    binCounts = hist(allTH,xBins);
    maxCounts = max(binCounts);
    bar(xBins,binCounts);
    hold on
    hl = line([hdi95(1),hdi95(2)],[0,0]);
    set(hl,'Color','y','LineWidth',10);
    axis([0,1,0,maxCounts]);
    ax = axis;
    xlabel('\theta');
    ylabel('Frequency');
    title('Sampled distribution');
    tStr = sprintf('Proposal SD = %0.2f\nMode = %0.3f\n95%%HDI = %0.3f,%0.3f',...
        mySigma,myMode,hdi95(1),hdi95(2));
    text(0.1,0.6*(ax(4)-ax(3)),tStr);
    
    % superimpose the correct answer
    y = betapdf(xBins,a+z,b+N-z);
    % need to scale for histogram
    y = (y ./ max(y)) .* ax(4);
    plot(xBins,y,'r-');
    
    subplot(3,1,2);
    plot(allTH(end-100:end),nSims-100:nSims,'bo-');
    hold on
    ax = axis;
    axis([0,1,ax(3),ax(4)]);
    xlabel('\theta');
    ylabel('Step in chain');
    title('End of chain');
    tStr = sprintf('N_a_c_c / N_p_r_o = %0.3f',(nAccepted/nSims));
    text(0.1,ax(3) + (0.6*(ax(4)-ax(3))),tStr);
    
    subplot(3,1,3);
    plot(allTH(1:100),1:100,'bo-');
    ax = axis;
    axis([0,1,ax(3),ax(4)]);
    xlabel('\theta');
    ylabel('Step in chain');
    title('Beginning of chain');
    

end

