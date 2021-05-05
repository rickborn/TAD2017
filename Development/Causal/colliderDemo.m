% colliderDemo.m
%
% How conditioning on a collider introduces spurious correlations
% Inspired by Judea Pearl's "The Book of Why", see pp. 197-200.
%
% RTB wrote it, 12 January 2021; post Trump mob storming capital

%% Generate some fake data for GPA and GRE scores that are uncorrelated:

nApps = 100;    % # of applications to grad school

% GPAs, uniformly distributed over [2,4):
gpa = (rand(nApps,1) .* 2) + 2;

% GRE scores normally distributed around a mean of 600, with s.d. of 40:
gre = normrnd(600,40,nApps,1);

figure
plot(gre,gpa,'k+');
hold on
xlabel('GRE Score (quant)'); ylabel('GPA (science)');
% plot least squares regression line for each data set
lsline

% calculate correlation:
[rhoAll,pvalAll] = corr(gre,gpa);

%% Decide who gets accepted: 
% normalize each metric, add them together, then accept the top half

nlGPA = (gpa - min(gpa)) ./ (max(gpa) - min(gpa));
nlGRE = (gre - min(gre)) ./ (max(gre) - min(gre));
totScore = nlGPA + nlGRE;
acc = totScore > median(totScore);

%% Now look to see if we have introduced a correlation:

plot(gre(acc),gpa(acc),'ro');
%('GRE of those accepted'); ylabel('GPA of those accepted');
% plot least squares regression line for each data set
lsline

% calculate correlation:
[rhoAcc,pvalAcc] = corr(gre(acc),gpa(acc));
tStr = sprintf('R_a_l_l=%0.2f, P_a_l_l= %0.3f, R_a_c_c=%0.2f, P_a_c_c=%0.3f',...
    rhoAll,pvalAll,rhoAcc,pvalAcc);
title(tStr);

legend('All','All','Accepted','Accepted','Location','NorthEast');

% %% Use "real" data, where we start with pos. corr.
% 
% ds82 = readtable('Grad_School_82.xlsx'); % All graduate programs (*census*)
% ds15 = readtable('Grad_School_15.xlsx'); % random sample of 15
% 
% %% Plot GPA and GRE scores
% 
% figure
% plot(ds82.GRE,ds82.GPA,'k+');
% hold on
% plot(ds15.GRE,ds15.GPA,'ro');
% xlabel('GRE Score (quant)'); ylabel('GPA (science)');
% % plot least squares regression line for each data set
% lsline
% legend('Census','Sample','Sample','Census','Location','NorthWest');
% 
% [rho,pval] = corr(ds82.GRE,ds82.GPA);
% tStr = sprintf('R = %0.2f, pVal = %0.3f',rho,pval);
% title(tStr);
% 
% %% Same drill: Decide who gets accepted: 
% % normalize each metric, add them together, then accept the top half
% 
% nlGPA = nl(ds82.GPA);
% nlGRE = nl(ds82.GRE);
% totScore = nlGPA + nlGRE;
% acc = totScore > median(totScore);
% 
% %% Condition on acceptance:
% 
% figure
% plot(ds82.GRE(acc),ds82.GPA(acc),'k+');
% xlabel('GRE of those accepted'); ylabel('GPA of those accepted');
% % plot least squares regression line for each data set
% lsline
% 
% [rho,pval] = corr(ds82.GRE(acc),ds82.GPA(acc));
% tStr = sprintf('R = %0.2f, pVal = %0.3f',rho,pval);
% title(tStr);
