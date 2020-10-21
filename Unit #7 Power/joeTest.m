% joeTest.m
%
% Testing Joe Hurley's code

rng default;
alpha = 0.05;
n_sims = 10000;
n_start = 20;
n_increment = 5;
n_end = 50;
group_expt = randn(n_end, n_sims);
group_ctrl = randn(n_end, n_sims);
sig_results = zeros(n_sims, 1);
p_results = zeros(n_sims, 1);
n_results = zeros(n_sims, 1);

for col = 1:n_sims
  sig_reached = boolean(0);
  n_current = n_start;
  while sig_reached == boolean(0) && n_current <= n_end
    [~, p] = ttest2(group_expt(1:n_current, col), group_ctrl(1:n_current, col));
    if p <= alpha
      sig_reached = boolean(1);
    end
    n_current = n_current + n_increment;
  end
  sig_results(col) = sig_reached;
  p_results(col) = p;
  n_results(col) = n_current - n_increment;
end
false_pos_rate = sum(sig_results) / n_sims * 100;

% So, with n_sims = 10,000, the result is 12.5%, but
% with n_sims = 100,000, the result is 13.24%. It all has to do with the
% order in which the random arrays are filled in! Interesting lesson here.