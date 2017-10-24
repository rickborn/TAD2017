% FWERwithCorrSim.m
%
% solution code for Thornquist test question

myCorr = 0.6;   % correlation among tests
nSamp = 20;     % number per group
muExpt = [1,1,1,1,1];
muCtrl = [0,0,0,0,0];
mySigma = (ones(5,5) .* myCorr) + (eye(5) .* (1 - myCorr));

% What is the probability of at least one false positive?
nSim = 10000;
nFP = 0;
for k = 1:nSim
    R1 = mvnrnd(muCtrl,mySigma,nSamp);
    R2 = mvnrnd(muCtrl,mySigma,nSamp);
    H = ttest2(R1,R2);
    if any(H)
        nFP = nFP + 1;
    end
end
pFP = nFP / nSim;

% What is the probability of producing a false positive if the tests are
% totally uncorrelated (i.e. sigma is 0 everywhere except the diagonal)?
mySigma = eye(5);   % and repeat the above

% RTB Question: What is the probability if all tests are perfectly
% correlated? (intuition here is that you are now really only doing one
% test, so it should be close to 0.05.
mySigma = ones(5,5);
% Yup! I got 0.0529
