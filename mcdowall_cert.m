%% Sovereign Default and Capital Controls
% Matlab File 1: certainty case
% By: Robert A. McDowall

clear all
clc
global y_0 y_1 g_0 gamma R tau phi; 

y_0 = 1.05;
y_1 = 1;
g_0 = .35;
beta = .96;
R = 1.04;
q = 1/R; 
gamma = 2;
phi = 0.12; 
%%

util = @(c) (c.^(1-gamma))./(1-gamma);
mutil = @(c) c.^(-gamma);


%% Comparison - Govt. Spending

g_vec = 0:.01:.6;
iter = 5;
B_trans = 0;


for i=1:length(g_vec)
    
    g_0 = g_vec(i);
       
    funsolve = @(B_d) -(util(y_0 -q*B_d) +beta*util(y_1 +B_d - (g_vec(i)*R)));
    [Bd_max] = fmincon(funsolve,.05, [], [], [], [], 0 ,g_0/R);

    Bd_vec(i) = Bd_max;
    B_vec(i) = g_0*R;
    
    Bf_vec(i) = B_vec(i) - Bd_max;
    
    c_0 = y_0 - q*Bd_max;
    c_1 = y_1 + Bd_max-B_vec(i);
    optc_vec(i) = beta*mutil(c_1) - q* mutil(c_0);
    
    util_vec(i) = util(c_0) + beta* util(c_1);
    
end

%%
diff = util(y_1 + Bd_vec-B_vec) - util(y_1 - phi);
%%
close all
plot(g_vec, util_vec)

%% Find Autarchy Values

for i = 1:length(diff)
    
    if diff(i) > 0
        
        util_final(i) = util_vec(i);
        
    else
        
        util_final(i) =  util(y_0 - g_vec(i)) + beta* util(y_1);
        
    end
    
    
end
%%
close all
plot(g_vec(21:end), util_final(21:end))
hold on 
plot(g_vec, util_vec, 'green')
hold on 
plot(g_vec(1:20), util_final(1:20))

%% Opitmal Controls
for i = 1:length(diff)
    
    if diff(i) > 0
        
        util_CC(i) = util_vec(i);
        
    else
      
      %Solve optimal savings implied by optimal tau
      tau_s = 1 - 1/(( mutil(y_0 - g_vec(i) + phi/R))/(beta*R*mutil(y_1 - phi)));
      funsolve = @(B_d) -(util(y_0 - q*B_d*(1-tau_s)) +beta*util(y_1 + B_d - (g_vec(i)*R+tau_s*B_d)));
      
      [Bd_max] = fmincon(funsolve,.1, [], [], [], [], 0 ,g_vec(i)/R);
      
      Bd_CC(i) = Bd_max; 
      B_CC(i) = (g_vec(i)*R + tau_s*Bd_CC(i));
       
      c_0 = y_0 - q*Bd_CC(i)*(1-tau_s);
      c_1 = y_1 + Bd_CC(i)-B_CC(i);
        
      util_CC(i) =  util(c_0) + beta* util(c_1);
    end
       
end

%% Plot
repay_index = max(find(diff>0))
 
close all
plot(g_vec/y_0, util_vec, 'green', 'LineWidth',2.1)
hold on
plot(g_vec(repay_index:end)/y_0, util_final(repay_index:end), '*b', 'LineWidth',1)
axis([.15,.58,-3.2,-1.9])
hold on 
plot(g_vec/y_0, util_CC, '--red', 'LineWidth',2.2) 
hold on
plot(g_vec(1:repay_index)/y_0, util_vec(1:repay_index), '*b', 'LineWidth',1 )
title('Optimal Capital Controls','interpreter','latex','fontsize',18)
ylabel('Utility','interpreter','latex','fontsize',18)
xlabel('Government Spending (\% GDP)','interpreter','latex','fontsize',18)
legend('Commitment', 'No Commitment', 'Capital Controls')% ,14,'Location','Northeast')

%% Summary Statistics
clc 
index = 36;
g_val = g_vec(index)

%Internal debts
first_best_prop =Bd_vec(index)/B_vec(index)
ccont_prop = Bd_CC(index)/B_CC(index)

%convert from subsidy to tax space
tau_tbl = ( mutil(y_0 - g_val + phi/R))/(beta*R*mutil(y_1 - phi)) - 1 

%Welfare
util_vec(index) - util_CC(index)
util_vec(index) / util_CC(index)
