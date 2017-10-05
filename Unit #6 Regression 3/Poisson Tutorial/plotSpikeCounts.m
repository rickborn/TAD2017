function plotSpikeCounts(counts)

t = counts';

%figure(3);
subplot(2, 1, 1);
hist(t, 0:max(max(t)) * 1.1);
s = axis;
axis([0 s(2) s(3) s(4)]);
xlabel('Spike Counts');
ylabel('#');

subplot(2, 1, 2);
m = mean(t);
v = var(t);
plot(m, v, 'bo');
s = axis;
axisLimit = max(s(2), s(4));
axis([0 axisLimit 0 axisLimit]);
hold on;
line([0, axisLimit], [0, axisLimit]);
hold off;
axis('square');
xlabel('Mean Spike Count');
ylabel('Spike Count Variance');

end