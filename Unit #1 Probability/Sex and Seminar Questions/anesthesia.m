% anesthesia.m
%
% In an experiment with anesthesia, subjects were asked to "respond" under
% different anesthetic regimes, tracked with EEG over the frontal lobes and
% analyzed by the phase of alpha waves w/r/t slow delta waves, stratified
% into 2 groups: 'troughMax' meaning the alpha was concentrated in the
% troughs of the delta waves (postulated to indicate lighter anesthesia)
% and 'peakMax'.
%
% RTB wrote it for Emery Brown's SAC meeting, 03 August 2018

% underlying probability of a subject responding
allP = 0:0.01:1;
% flat prior
prior = ones(size(allP)) ./ length(allP);

% in the 'peakMax' group:
likePM = binopdf(0,50,allP);
% and the 'troughMax' group:
likeTM = binopdf(6,36,allP);

% calculate the posterior probabilities
postPM = likePM .* prior;
postPM = postPM ./ trapz(postPM);
postTM = likeTM .* prior;
postTM = postTM ./ trapz(postTM);

% plot 'em
figure, plot(allP,postTM,'k-')
hold on
plot(allP,postPM,'r-')
legend('troughMax','peakMax')
xlabel('Probability of response');
ylabel('Posterior probability')