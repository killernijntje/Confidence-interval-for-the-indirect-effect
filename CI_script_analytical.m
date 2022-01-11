clc;clear;
tic 
%SEED
rng(123);
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
Beta_true  = 1; 
Tau_true   = 1; 
Gamma_true = Alpha_true * Beta_true;
Chi_crit= 2.064^2;
scale_param=1;



    % DGP sample emulating

    X = unifrnd(0,1,N,1);
    X = X - mean(X);
    esp1= normrnd(0,1,N,1);
    esp2= normrnd(0,1,N,1);
    
    M = Alpha_true .* X  + esp2;
    Y = Tau_true .* X + Beta_true.*M + esp1;
    
    XM=[X M];
    
    % estimating Tau, Aplha and Beta 
   
    Alpha_hat = inv(X'*X).*X'*M;
    Equation1_hats = (XM'*XM)\XM'*Y; 
    Beta_hat = Equation1_hats(2);
    Tau_hat = Equation1_hats(1);
    
    %estimating Y and M 
    
    Y_hat = Tau_hat.*X + Beta_hat.*M; 
    M_hat = Alpha_hat.*X;
    
    % estimating variance of ESP 1 and ESP 2

    esp1_variance_hat = (Y-Y_hat)'*(Y-Y_hat)/(N-k1);
    esp2_variance_hat = (M-M_hat)'*(M-M_hat)/(N-k2);

    % Estimating variance of Alpha_hat and Beta_hat 
    
    Beta_hat_variance  = inv(M'*M).*esp1_variance_hat;
    Alpha_hat_variance = inv(X'*X).*esp2_variance_hat;
 
   %Calculating the bounds for Alpha
    
    Alphas_sol = zeros(1,4);
    Betas_sol = zeros(1,4);
    
    Alphas_sol(1)= a1(Alpha_hat, Beta_hat, sqrt(Alpha_hat_variance), sqrt(Beta_hat_variance), scale_param, Chi_crit);
    Alphas_sol(2)= a2(Alpha_hat, Beta_hat, sqrt(Alpha_hat_variance), sqrt(Beta_hat_variance), scale_param, Chi_crit);
    Alphas_sol(3)= a3(Alpha_hat, Beta_hat, sqrt(Alpha_hat_variance), sqrt(Beta_hat_variance), scale_param, Chi_crit);
    Alphas_sol(4)= a4(Alpha_hat, Beta_hat, sqrt(Alpha_hat_variance), sqrt(Beta_hat_variance), scale_param, Chi_crit);

    %Calculating the bounds for Beta
    
    Betas_sol(1)= b1(Alpha_hat, Beta_hat, sqrt(Alpha_hat_variance), sqrt(Beta_hat_variance), scale_param, Chi_crit);
    Betas_sol(2)= b2(Alpha_hat, Beta_hat, sqrt(Alpha_hat_variance), sqrt(Beta_hat_variance), scale_param, Chi_crit);
    Betas_sol(3)= b3(Alpha_hat, Beta_hat, sqrt(Alpha_hat_variance), sqrt(Beta_hat_variance), scale_param, Chi_crit);
    Betas_sol(4)= b4(Alpha_hat, Beta_hat, sqrt(Alpha_hat_variance), sqrt(Beta_hat_variance), scale_param, Chi_crit);

    %Calculating Gamma values 
    %removing complex values with img part <10^-8
    
    Gammas_sol = Alphas_sol.*Betas_sol ;
    Gammas_sol = real(Gammas_sol(abs(imag(Gammas_sol)) < 10^(-8)));
    CI_est = [min(Gammas_sol) max(Gammas_sol)] ;
    
    


    
%     
%     
%     
%         %calculating all functions needed for Alpha
%     a_func_1 = a_func_1(Alpha_hat, Beta_hat, sqrt(Alpha_hat_variance), sqrt(Beta_hat_variance), scale_param, Chi_crit);
%     a_func_2_inner = a_func_2_inner(Alpha_hat, Beta_hat, sqrt(Alpha_hat_variance), sqrt(Beta_hat_variance), scale_param,Chi_crit);
%     a_func_2_p = a_func_2_p(Alpha_hat, Beta_hat, sqrt(Alpha_hat_variance), sqrt(Beta_hat_variance), scale_param,Chi_crit,a_func_2_inner);
%     a_func_2_m = a_func_2_m(Alpha_hat, Beta_hat, sqrt(Alpha_hat_variance), sqrt(Beta_hat_variance), scale_param,Chi_crit,a_func_2_inner);
%     
%     
%     %Calculating all 4 solutions for Alpha
%     
%     alphas_sol(1)=3/4 * Alpha_hat - a_func_1 - a_func_2_m ; 
%     alphas_sol(2)=3/4 * Alpha_hat - a_func_1 + a_func_2_m ; 
%     alphas_sol(3)=3/4 * Alpha_hat + a_func_1 - a_func_2_p ; 
%     alphas_sol(4)=3/4 * Alpha_hat + a_func_1 + a_func_2_p ; 
%     
%     
    
    
%     b_func_inner_1 = b_func_1_inner_1(Alpha_hat, Beta_hat, sqrt(Alpha_hat_variance), sqrt(Beta_hat_variance), scale_param, Chi_crit);
%     b_func_inner_2 = b_func_1_inner_2(Alpha_hat, Beta_hat, sqrt(Alpha_hat_variance), sqrt(Beta_hat_variance), scale_param, Chi_crit);
%     b_func_p = b_func_1(Alpha_hat, Beta_hat, sqrt(Alpha_hat_variance), sqrt(Beta_hat_variance), scale_param, Chi_crit,b_func_inner_1,b_func_inner_2);
%     b_func_m = b_func_1(Alpha_hat, Beta_hat, sqrt(Alpha_hat_variance), sqrt(Beta_hat_variance), scale_param, Chi_crit,-1*b_func_inner_1,-1*b_func_inner_2);
%     
%     
%     betas_sol(1)= ( 3/4 * Beta_hat^2 * Alpha_hat_variance + b_func_m)  / (Beta_hat*Alpha_hat_variance) ; %correct
%     be        tas_sol(2)= ( 3/4 * Beta_hat^2 * Alpha_hat_variance - b_func_m)  / (Beta_hat*Alpha_hat_variance) ;
%     betas_sol(3)= ( 3/4 * Beta_hat^2 * Alpha_hat_variance - b_func_p)  / (Beta_hat*Alpha_hat_variance) ;
%     betas_sol(4)= ( 3/4 * Beta_hat^2 * Alpha_hat_variance + b_func_p)  / (Beta_hat*Alpha_hat_variance) ; %correct
%  
   