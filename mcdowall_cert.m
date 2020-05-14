clear all
clc
global y_0 y_1 g_0 gamma R tau phi; 

y_0 = 1.05;
y_1 = 1;
g_0 = .35;
beta = .96;
R = 1.04;
q = 1/R; 
tau = 0;
gamma = 2;
phi = 0.12; 
%%

util = @(c) (c.^(1-gamma))./(1-gamma);
mutil = @(c) c.^(-gamma);

%util = @(c) log(c);
%mutil = @(c) 1./c;

funsolve = @(B_d) -(util(y_0 -q*B_d*(1+tau)) +beta*util(y_1 +B_d - (g_0-q*tau*B_d)/q));


%%
[Bd_max] = fmincon(funsolve,.05, [], [], [], [], 0 ,g_0/R)

%% Comparison - Govt. Spending

g_vec = 0:.01:.6;
iter = 5;
B_trans = 0;
tau = 0; 

for i=1:length(g_vec)
    
    g_0 = g_vec(i);
    
    for j = 1:iter
        
    funsolve = @(B_d) -(util(y_0 -q*B_d*(1+tau)) +beta*util(y_1 +B_d - (g_vec(i)-q*tau*B_trans)/q));
    [Bd_max] = fmincon(funsolve,.05, [], [], [], [], 0 ,g_0/R);
    
    B_trans = Bd_max;

    end
    
    Bd_vec(i) = Bd_max;
    B_vec(i) = (g_0-q*tau*Bd_max)/q;
    
    Bf_vec(i) = B_vec(i) - Bd_max;
    
    c_0 = y_0 - q*Bd_max*(1+tau);
    c_1 = y_1 + Bd_max-B_vec(i);
    optc_vec(i) = beta*mutil(c_1) - q* mutil(c_0)*(1+tau);
    
    util_vec(i) = util(c_0) + beta* util(c_1);
    
end

%%

diff = util(y_1 + Bd_vec-B_vec) - util(y_1 - phi);
%%
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
plot(g_vec(21:end), util_final(21:end))
hold on 
plot(g_vec, util_vec, 'green')
hold on 
plot(g_vec(1:20), util_final(1:20))

%% Opitmal Controls

for i = 1:length(diff)
    
    if diff(i) > 0
        
        util_CC(i) = util_vec(i)
        
    else
        
        %search for fixed point on tau that satisfied IC constraint
        tau = -0.1; 
        g_0 = g_vec(i);
        B_trans = 0.1;
       
       for k = 1:iter 
            for j = 1:iter
        
                funsolve = @(B_d) -(util(y_0 -q*B_d*(1+tau)) +beta*util(y_1 +B_d - (g_vec(i)-q*tau*B_trans)/q));
                [Bd_max] = fmincon(funsolve,.1, [], [], [], [], 0 ,g_0/R);
    
                 B_trans = Bd_max;

            end
        
           Bd_CC(i) = Bd_max;
            tau_opt = (g_0/q - phi)/Bd_max-1; 
        
            tau = tau_opt; 
       end 
        
       tau_optvec(i) = tau_opt;  
       
       B_CC(i) = (g_0-q*tau*Bd_max)/q;
       
       c_0 = y_0 - q*Bd_max*(1+tau_opt);
       c_1 = y_1 + Bd_max-B_CC(i);
        
       util_CC(i) =  util(c_0) + beta* util(c_1);
    end
    
    
end
 %%
 close all
plot(g_vec/y_0, util_vec, 'green', 'LineWidth',2.1)
hold on
plot(g_vec(41:end)/y_0, util_final(41:end), '*b', 'LineWidth',1)
axis([.15,.58,-3.2,-1.9])
hold on 
plot(g_vec/y_0, util_CC, '--red', 'LineWidth',2.2) 
hold on
plot(g_vec(1:41)/y_0, util_vec(1:41), '*b', 'LineWidth',1 )
title('Optimal Capital Controls','interpreter','latex','fontsize',18)
ylabel('Utility','interpreter','latex','fontsize',18)
xlabel('Government Spending (\% GDP)','interpreter','latex','fontsize',18)
legend('Commitment', 'No Commitment', 'Capital Controls' ,14,'Location','Northeast')

%% Summary Statistics

Bd_CC(26)/B_CC(26)
Bd_vec(26)/B_vec(26)

tau_optvec(26)




