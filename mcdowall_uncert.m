%% Sovereign Debt and Capital Controls
% Matlab File 2: uncertainty case
% By: Robert A. McDowall
clear all
clc
global y_0 y_1 g_0 beta R;  
%% Parameters and Functional Forms
y_0 =1.05;
y_1 = 1;
g_0 = .3;
R = 1.04;

beta = .96;

util = @(c) log(c);
mutil = @(c) 1./c;


%%Log Normal Distribution, mean preserving spread
mu = .175;
sigma = [.001:.0001:.08];

mu_vec = log((mu^2)./sqrt(sigma + mu.^2));

sigma_vec = sqrt(log(sigma./mu^2 + 1));


%% Solve on grid

for j = 1:length(sigma_vec)
  
sigma = sigma_vec(j);
mu = mu_vec(j);
    
tau_veca = 0.001:.0002:.3;


for i = 1:length(tau_veca)
         
    tau = tau_veca(i); 
    
   solveBf = @(Bf)  -Bf*(1+tau) + g_0/((1-logncdf(Bf,mu,sigma))/(R * (1+tau))) ...
    - ((beta* (1-logncdf(Bf,mu,sigma))*y_0)...
    - (((1-logncdf(Bf,mu,sigma))/(R * (1+tau)))*(y_1 + tau*Bf) - g_0))/...
    (((1-logncdf(Bf,mu,sigma))/(R * (1+tau)))*(beta* (1-logncdf(Bf,mu,sigma)) + 1));

    options = optimset(optimset('fsolve'), 'TolFun', 1.0e-6, 'TolX', 1.0e-6);
    options = optimset('Display','off');

    Bf_solve = fsolve(solveBf, .1, options);  
    
    Bf_vec2(i,j) = max(0,Bf_solve);
    
end
 
end

%% Calculate Solutions
tau_vec2 = repmat(tau_veca,length(sigma_vec),1)';


for i = 1:length(sigma_vec)
    
    F_vec(:,i) = max(0, logncdf(Bf_vec2(:,i),mu_vec(i),sigma_vec(i))); 

end
  
q_vec2 = (1-F_vec)/(R*(1+tau));

Bd_vec2 = (g_0./q_vec2) - Bf_vec2.*(1+tau_vec2);

Rev = q_vec2.*(Bf_vec2 + Bd_vec2) + q_vec2.*(tau_vec2).* Bf_vec2;

B_ratio = Bd_vec2 ./ (Bd_vec2+Bf_vec2);

%% Find Optima

c_0V = y_0 - q_vec2.*Bd_vec2;

c_1rV = y_1 - Bf_vec2;

expect_dV = @(phi) util(y_1 - phi)*lognpdf(mu, mu, sigma);

for i = 1:length(sigma_vec)
    
    sigma = sigma_vec(i);
    mu = mu_vec(i);
    expect_dV = @(phi) util(y_1 - phi)*lognpdf(mu, mu, sigma);
    
    util_dV(i) = quadgk(expect_dV,0,y_1 - .2);

end

util_dV2 =repmat(util_dV, size(c_0V,1), 1);
%% Calculate Value Functions

util_m = [log(c_0V) + beta.*[(1-F_vec) .* log(c_1rV) + F_vec.* util_dV2]]

close all
%%
plot(tau_veca,util_m)

[M,I] = max(util_m);

var_phi = sigma_vec.^2;

%% Maxima (for plots)

for i=1:length(I)
    Rep_max(i) = 1-F_vec(I(i), i);
    q_max(i) = q_vec2(I(i), i);
    tau_max(i) = tau_vec2(I(i), i);
    B_ratmax(i) = B_ratio(I(i), i);
    Bf_max(i) = Bf_vec2(I(i), i);
    Bd_max(i) = Bd_vec2(I(i), i);
    
end

%%

figure

subplot(1,2,1)
[ax h1 h2] = plotyy(var_phi,tau_max, var_phi, B_ratmax);
title('A.) Optimal Controls', 'interpreter','latex','fontsize',16)
axes(ax(1)); ylabel('Optimal Control', 'interpreter','latex','fontsize',16);
axes(ax(2)); ylabel('Domestic Debt Share', 'interpreter','latex','fontsize',16);
xlabel('Variance of $\Phi$','interpreter','latex','fontsize',16);
set(ax,{'ycolor'},{'b';'r'})
%set(ax(1),'YLim',[-.52 -.36])
%set(ax(2),'YLim',[.76 .92])
%set(ax(2),'YTick',[0:.1:.55])
%set(ax(1),'YTick',[-.45:.01:-.4])
set(h1,'LineWidth', 2)
set(h2, 'LineStyle','--', 'color','red', 'LineWidth', 2)
set(ax,'xlim',[0,max(var_phi)]);

subplot(1,2,2)
[ax h1 h2] = plotyy(var_phi,q_max, var_phi, Bf_max);
title('B.) Bond Prices', 'interpreter','latex','fontsize',16)
axes(ax(1)); ylabel('Bond Prices, $q$', 'interpreter','latex','fontsize',16);
axes(ax(2)); ylabel('Foreign Lending, $B_f$', 'interpreter','latex','fontsize',16);
xlabel('Variance of $\Phi$','interpreter','latex','fontsize',16);
set(ax,{'ycolor'},{'k';'r'})
set(h1,'color','k', 'LineWidth', 2)
set(h2, 'LineStyle','--', 'color','r', 'LineWidth', 2)
set(ax,'xlim',[0,max(var_phi)]);
%set(ax(1),'YLim',[0 1.2])
%set(ax(2),'YLim',[.05 .12])

%%
BfDashed = Bf_max;
BfDashed(1:2:length(BfDashed)) = NaN;

var_phiDashed = var_phi;
var_phiDashed(1:4:length(var_phiDashed)) = NaN;

figure
title('Optimum', 'interpreter','latex','fontsize',14)

subplot(2,2,1)
plot(var_phi,tau_max, 'r--', 'LineWidth',2.25 )
title('Capital Controls', 'interpreter','latex','fontsize',14)
xlabel('Variance of $\Phi$','interpreter','latex','fontsize',12);
ylabel('$\tau$','interpreter','latex','fontsize',12);

subplot(2,2,2)
plot(var_phi, B_ratmax, 'b', 'LineWidth',2.25 )
title('Domestic Share', 'interpreter','latex','fontsize',14)
xlabel('Variance of $\Phi$','interpreter','latex','fontsize',12);
ylabel('$\frac{B_d}{B}$','interpreter','latex','fontsize',12);

subplot(2,2,3)
plot(var_phi,q_max, 'b--', 'LineWidth',2.25 )
title('Bond Prices', 'interpreter','latex','fontsize',14)
xlabel('Variance of $\Phi$','interpreter','latex','fontsize',12);
ylabel('$q$','interpreter','latex','fontsize',12);

subplot(2,2,4)
plot(var_phi, Bf_max,  'r', 'LineWidth',2.25)
title('Foreign Holdings', 'interpreter','latex','fontsize',14)
xlabel('Variance of $\Phi$','interpreter','latex','fontsize',12);
ylabel('$B_f$','interpreter','latex','fontsize',12);

%% Plots 2

figure

subplot(2,2,3)
plot(tau_vec2(:,5), q_vec2(:,5), 'LineWidth',2.25)
hold on
plot(tau_max(5),q_max(5),'bo','MarkerSize',10)
hold on
plot(tau_vec2(:,800), q_vec2(:,800), 'r--', 'LineWidth',2.25)
hold on
plot(tau_max(800), q_max(800),'ro','MarkerSize',10)
title('Bond Prices', 'interpreter','latex','fontsize',14)
ylabel('$q$', 'interpreter','latex','fontsize',14) % y-axis label
xlabel('$\tau$','interpreter','latex','fontsize',14);
xlim([0 .3])
%grid minor


subplot(2,2,1)
plot(tau_vec2(:,5), Bd_vec2(:,5), 'LineWidth',2.25)
hold on
plot(tau_max(5),Bd_max(5),'bo','MarkerSize',10)
hold on
plot(tau_vec2(:,800), Bd_vec2(:,800), 'r--', 'LineWidth',2.25)
hold on
plot(tau_max(800), Bd_max(800),'ro','MarkerSize',10)
%ylim([min(q_vec)-.03, max(q_vec)+.02])
title('Domestic Holdings', 'interpreter','latex','fontsize',14)
ylabel('$B_d$', 'interpreter','latex','fontsize',14) % y-axis label
xlabel('$\tau$','interpreter','latex','fontsize',14);
%grid minor
h_legend = legend('low var.', 'low opt.', 'high var.', 'high opt.');
set(h_legend,'interpreter','latex','FontSize',10);
xlim([0 .3])

subplot(2,2,4)
plot(tau_vec2(:,5), Bf_vec2(:,5), 'LineWidth',2.25)
hold on
plot(tau_max(5),Bf_max(5),'bo', 'MarkerSize',10)
hold on
plot(tau_vec2(:,800), Bf_vec2(:,800), 'r--', 'LineWidth',2.25)
hold on
plot(tau_max(800), Bf_max(800),'ro','MarkerSize',10)
title('Foreign Holdings', 'interpreter','latex','fontsize',14)
ylabel('$B_f$', 'interpreter','latex','fontsize',14) % y-axis label
xlabel('$\tau$','interpreter','latex','fontsize',14);
xlim([0 .3])
%grid minor

subplot(2,2,2)
plot(tau_vec2(:,5), B_ratio(:,5), 'LineWidth',2.25)
hold on
plot(tau_max(5),B_ratmax(5),'bo','MarkerSize',10)
hold on
plot(tau_vec2(:,800), B_ratio(:,800), 'r--', 'LineWidth',2.25)
hold on
plot(tau_max(800),B_ratmax(800),'ro', 'MarkerSize',10)
title('Domestic Share', 'interpreter','latex','fontsize',14)
ylabel('$\frac{B_d}{B}$', 'interpreter','latex','fontsize',14) % y-axis label
xlabel('$\tau$','interpreter','latex','fontsize',14);
xlim([0 .3])
%grid minor

%%
close all
plot(tau_max, Rep_max,  'r', 'LineWidth',2.25)
%title('Foreign Holdings', 'interpreter','latex','fontsize',14)
%xlabel('Variance of $\Phi$','interpreter','latex','fontsize',12);
%ylabel('$B_f$','interpreter','latex','fontsize',12);
%%
close all
%plot(tau_vec2(:,5), (1-F_vec(:,5)), 'LineWidth',2.25)
%hold on
%plot(tau_max(5),Rep_max(5),'ko','MarkerSize',10)
%hold on
%plot(tau_vec2(:,50), (1-F_vec(:,50)), 'r--', 'LineWidth',2.25)
%hold on
%plot(tau_max(50),Rep_max(50),'ko','MarkerSize',10)
%hold on
plot(tau_vec2(:,20), (1-F_vec(:,20)), 'r--', 'LineWidth',3)
hold on
plot(tau_max(20),Rep_max(20),'ko','MarkerSize',15, 'LineWidth',4)
title('Mitigation of Default Risk', 'interpreter','latex','fontsize',18)
ylabel('Repayment Probability','interpreter','latex','fontsize',16);
xlabel('$\tau$','interpreter','latex','fontsize',16);

%%
close all
plot(tau_vec2(:,990), Bf_vec2(:,990), 'r--', 'LineWidth',2.75)
hold on
plot(tau_max(990),Bf_max(990),'ko','MarkerSize',15)
title('Mitigation of Default Risk', 'interpreter','latex','fontsize',14)
ylabel('Repayment Probability','interpreter','latex','fontsize',12);
xlabel('$\tau$','interpreter','latex','fontsize',12);
