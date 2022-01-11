clc;

%alpha ranges from -1 to 1 
%beta ranges from 0 to 1 
%CI coverage is plotted in the Z-axis
%
%Scale parameter for the Chi critical value
%can be chosen by a number betwee   n 3 and 18
%
% 3=0.25    9=0.55      15=0.85
% 4=0.30    10=0.60     16=0.90
% 5=0.35    11=0.65     17=0.95
% 6=0.40    12=0.70     18=1
% 7=0.45    13=0.75
% 8=0.50    14=0.80

interpolated_values= zeros(1,5);
%%
for i=5:18 
scaleparameter=i;

%interpolating 

g=Scale_parameter_result(scaleparameter,2:232)';


ab=Scale_parameter_result(1:2,2:232)';


F.Method = 'natural';
F.ExtrapolationMethod = 'none'; 
F=scatteredInterpolant(ab,g); 
    
[xq,yq] = meshgrid(-0.5:0.001:0.5, 0:0.001:0.5);
vq = F(xq,yq);
    
%selecting values between 0.945 and 0.955

selector = (vq>0.948 & vq<0.952);

% changing matrices to collumn vectors
CI_interpolated_values_col= reshape(selector.*vq, [],1);
alpha_interpolated_values_col= reshape(selector.*(xq+10), [],1);
beta_interpolated_values_col= reshape(selector.*(yq +10), [],1);

% Elminitating zero entries
CI_interpolated_values_col(CI_interpolated_values_col==0)=[];
alpha_interpolated_values_col(alpha_interpolated_values_col==0)=[];
alpha_interpolated_values_col=alpha_interpolated_values_col-10;
beta_interpolated_values_col(beta_interpolated_values_col==0)=[];
beta_interpolated_values_col=beta_interpolated_values_col-10;
N_col = repmat(N,size(beta_interpolated_values_col)); 
detla_col = repmat(Scale_parameter_result(scaleparameter),size(beta_interpolated_values_col)); 

%creating data frame for regression

temp_inter_val = [CI_interpolated_values_col alpha_interpolated_values_col beta_interpolated_values_col N_col detla_col];
interpolated_values = [interpolated_values;temp_inter_val];
end 
%%

[add_a, add_b] = meshgrid(-1:0.1:1,0.5:0.1:1);
add_a = reshape(add_a,[],1);
add_b = reshape(add_b,[],1);
add_delta = ones(size(add_a));
add_CI = 0.95*add_delta; 
add_N = N * add_delta;

interpolated_values=[ interpolated_values; add_CI add_a add_b add_N add_delta  ] ;



%%

[add_a, add_b] = meshgrid(-10:0.1:-2,0:0.1:1);
add_a = reshape(add_a,[],1);
add_b = reshape(add_b,[],1);
add_delta = ones(size(add_a));
add_CI = 0.95*add_delta; 
add_N = N * add_delta;

interpolated_values=[ interpolated_values; add_CI add_a add_b add_N add_delta  ] ;

%%

[add_a, add_b] = meshgrid(2:0.1:10,0:0.1:1);
add_a = reshape(add_a,[],1);
add_b = reshape(add_b,[],1);
add_delta = ones(size(add_a));
add_CI = 0.95*add_delta; 
add_N = N * add_delta;

interpolated_values=[ interpolated_values; add_CI add_a add_b add_N add_delta  ] ;
%%

writematrix(interpolated_values,'interpolated_values_t.csv','Delimiter', ';');

%%
% 
% scaleparameter=i;
% %interpolating 
% 
% g=Scale_parameter_result(scaleparameter,2:232)';
% 
% 
% ab=Scale_parameter_result(1:2,2:232)';
% 
% 
% F.Method = 'natural';
% F.ExtrapolationMethod = 'none'; 
% F=scatteredInterpolant(ab,g); 
%     
% [xq,yq] = meshgrid(-1:0.001:1, 0:0.001:1);
% vq = F(xq,yq);
% 
% %selecting values between 0.945 and 0.955
% 
% selector = (vq>0.948 & vq<0.952);
% 
% % changing matrices to collumn vectors
% CI_interpolated_values_col= reshape(selector.*vq, [],1);
% alpha_interpolated_values_col= reshape(selector.*(xq+10), [],1);
% beta_interpolated_values_col= reshape(selector.*(yq +10), [],1);
% 
% % Elminitating zero entries
% CI_interpolated_values_col(CI_interpolated_values_col==0)=[];
% alpha_interpolated_values_col(alpha_interpolated_values_col==0)=[];
% alpha_interpolated_values_col=alpha_interpolated_values_col-10;
% beta_interpolated_values_col(beta_interpolated_values_col==0)=[];
% beta_interpolated_values_col=beta_interpolated_values_col-10;
% N_col = repmat(N,size(beta_interpolated_values_col)); 
% detla_col = repmat(Scale_parameter_result(scaleparameter),size(beta_interpolated_values_col)); 
% 
% %creating data frame for regression
% 
% temp_inter_val = [CI_interpolated_values_col alpha_interpolated_values_col beta_interpolated_values_col N_col detla_col];
% interpolated_values = [interpolated_values;temp_inter_val];
% 


%%
clc;
% fitting curve
fun1=@(c,interpolated_values)-0.6*(c(2)).^(c(3)*interpolated_values(:,2).^2+c(4)*interpolated_values(:,3).^2) +1 ;
fun2=@(c,interpolated_values)-0.6*exp(c(3)*interpolated_values(:,2).^2+c(4)*interpolated_values(:,3).^2) +1 ;

fun1_n=@(c,interpolated_values)-0.6*(c(2)).^(c(3)*interpolated_values(:,2).^2+c(4)*interpolated_values(:,3).^2 + c(5)*interpolated_values(:,2).^2.*interpolated_values(:,4)+c(6)*interpolated_values(:,3).^2.*interpolated_values(:,2)) +1 ;
fun2_n=@(c,interpolated_values)-0.6*exp(c(3)*interpolated_values(:,2).^2+c(4)*interpolated_values(:,3).^2 + c(5)*interpolated_values(:,2).^2.*interpolated_values(:,4)+c(6)*interpolated_values(:,3).^2.*interpolated_values(:,2)) +1 ;

% [c1,resnorm1] = lsqcurvefit(fun2, [ -1 1 -1 -1 -1] ,interpolated_values,interpolated_values(:,5)); 
% [c2,resnorm2] = lsqcurvefit(fun6, [ -1 1 -1 -1 -1] ,interpolated_values,interpolated_values(:,5)); 
[c1_n,resnorm1_n] = lsqcurvefit(fun1_n, [ -1 1 -1 -1 -1 -1] ,interpolated_values,interpolated_values(:,5)); 
[c2_n,resnorm2_n] = lsqcurvefit(fun2_n, [ -1 1 -1 -1 -1 -1] ,interpolated_values,interpolated_values(:,5)); 

