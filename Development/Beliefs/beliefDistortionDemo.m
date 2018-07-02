% beliefDistortionDemo.m
%
% reproduces fig. 2 of Zhang & Maloney 2012

trueP = 0.01:0.01:0.99;

figure;
subplot(1,2,1);
line([0,1],[0,1],'Color','k');
xlabel('p');
ylabel('\pi');
title('varying \gamma');
hold on;

allGamma = [0.2:0.4:1.8];
for k = 1:length(allGamma)
    q = [0.4,allGamma(k)];
    subjP = llo(trueP,q);
    plot(trueP,subjP,'b-');
end

subplot(1,2,2);
line([0,1],[0,1],'Color','k');
xlabel('p');
ylabel('\pi');
title('varying p_0');
hold on;

allP0 = [0.1:0.2:0.9];
for k = 1:length(allP0)
    q = [allP0(k),0.6];
    subjP = llo(trueP,q);
    plot(trueP,subjP,'b-');
end