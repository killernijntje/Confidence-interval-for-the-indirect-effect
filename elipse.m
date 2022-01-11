function [c,ceq] = elipse(x, Alpha_hat, Alpha_hat_variance, Beta_hat, Beta_hat_variance, C_sq_crit, N)
c = (x(1) -Alpha_hat)^2./(Alpha_hat_variance) + (x(2)-Beta_hat)^2./Beta_hat_variance - (-0.6*exp(-0.0157.*x(1).*N -0.21860.*x(2).*N)+1)*C_sq_crit;
ceq = [ ];