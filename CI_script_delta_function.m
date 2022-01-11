clc;clear;
tic 
%SEED
%rng(4)

%Simulating
%
% Y = tau * X + beta*M + eps1
% M = alpha * X + eps2
%
% Assumption: (eps1,eps2)' ~ N ( 0, diag(1,1))
%
% X is drawn UNIF(0,1) since X is exogenous 

% True parameters for Tau, Beta and Alpha can be changed. 
N  = 25;
k1 = 2; 
k2 = 1;
Alpha_true = 1;
Beta_true  = 0; 
Tau_true   = 1; 
Gamma_true = Alpha_true * Beta_true;

Sign_level= 0.05;
Crit_T = CritT(Sign_level, min([100 N-k1]), 'two'); %Calculates T_(alpha/2, N-K)
C_sq_crit = 3.84;
Z_crit = 1.96; % used to determine starting point x0 and plot axis
Chi_crit = C_sq_crit;
options = optimoptions(@fmincon, 'Display' , 'off');

    % DGP sample emulating

    X = unifrnd(0,1,N,1);
    X = X - mean(X);
    esp1= normrnd(0,1,N,1);
    esp2= normrnd(0,1,N,1);

    M = Alpha_true .* X  + esp1;
    Y = Tau_true .* X + Beta_true.*M + esp2;

    
        XM=[X M];
    
    % estimating Tau, Aplha and Beta 
   
    Alpha_hat = inv(X'*X).*X'*M;
    Equation1_hats = (XM'*XM)\XM'*Y; 
    Beta_hat = Equation1_hats(2);
    Tau_hat = Equation1_hats(1);
   
    Y_hat = Tau_hat.*X + Beta_hat.*M; 
    M_hat = Alpha_hat.*X;

    % estimating variance of ESP 1 and ESP 2

    esp1_variance_hat = (Y-Y_hat)'*(Y-Y_hat)/(N-k1);
    esp2_variance_hat = (M-M_hat)'*(M-M_hat)/(N-k2);

    % Estimating variance of Alpha_hat and Beta_hat 

    Alpha_hat_variance = inv(X'*X).*esp2_variance_hat;
    Beta_hat_variance  = inv(M'*M).*esp1_variance_hat;

    %declaring the optimization variables
    %a=alpha, b=beta

    a = optimvar('a');
    b = optimvar('b');

    % starting values 

    x0.a = [Alpha_hat+Z_crit*Alpha_hat_variance;Alpha_hat-Z_crit*Alpha_hat_variance];
    x0.b = [Beta_hat+Z_crit*Beta_hat_variance;Beta_hat-Z_crit*Beta_hat_variance];

    % defining the objective funtion and contraint
    gamma_obj = fcn2optimexpr(@objfunx,a,b,'OutputSize',[1,1]);
    elipse_const = (a -Alpha_hat)^2./(Alpha_hat_variance) + (b-Beta_hat)^2./Beta_hat_variance - C_sq_crit == 0 ; 
    
    %minimization for gamma_low
    gamma_low = optimproblem('Objective', gamma_obj);
    gamma_low.Constraints.elipse = elipse_const; 

    gamma_low_sol=struct('a', [0 , 0] ,'b', [0 , 0]);
    gamma_low_val=zeros(length(x0.b),length(x0.a));
    
    for j = 1:length(x0.b)
        for i = 1:length(x0.a)
            x0_temp.a=x0.a(i);
            x0_temp.b=x0.b(j);
           [gamma_low_sol_temp,gamma_low_val_temp] =  solve(gamma_low,x0_temp,'Options',options);
           gamma_low_sol(i,j)=gamma_low_sol_temp;
           gamma_low_val(i,j)=gamma_low_val_temp;
        end
    end 

    %maximization for gamma_high
    gamma_high = optimproblem('Objective', gamma_obj, 'ObjectiveSense','maximize');
    gamma_high.Constraints.elipse = elipse_const; 
    
    gamma_high_sol=struct('a', [0 , 0] ,'b', [0 , 0]);
    gamma_high_val=zeros(length(x0.b),length(x0.a));
    
    for j = 1:length(x0.b)
        for i = 1:length(x0.a)
            x0_temp.a=x0.a(i);
            x0_temp.b=x0.b(j);
           [gamma_high_sol_temp,gamma_high_val_temp] =  solve(gamma_high,x0_temp,'Options',options);
           gamma_high_sol(i,j)=gamma_high_sol_temp;
           gamma_high_val(i,j)=gamma_high_val_temp;
        end
    end

    %determining the lower and higher bound for Gamma (1-a)%CI
    [gamma_highest, gam_high_index]=max(gamma_high_val(:));
    [gamma_lowest, gam_low_index]=min(gamma_low_val(:));
    CI_est = [gamma_lowest, gamma_highest];



   
    
    
    
%creating plot of optimization problem
f = @objfunx;
g = @(a,b) (a -Alpha_hat)^2./(Alpha_hat_variance) + (b-Beta_hat)^2./Beta_hat_variance- C_sq_crit;
rnge = [-5 5 -5 5];
fimplicit(g,'k-'); 
axis(rnge);
hold on
fcontour(f, 'levellist', [gamma_lowest, gamma_highest]);
legend('Bivariate elipse','Gamma const.');
hold off


toc
