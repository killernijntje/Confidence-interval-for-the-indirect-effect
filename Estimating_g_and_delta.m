clc;clear;

fun1=@(c,interpolated_values)-0.6*exp(c(1)*interpolated_values(:,2).^2+c(2)*interpolated_values(:,3).^2) +1 ;
fun2=@(c,interpolated_values) -0.6*exp(c(1)*interpolated_values(:,2).^2.*interpolated_values(:,4)+c(2)*interpolated_values(:,3).^2.*interpolated_values(:,4)) +1 ;

fun2_lowerbound_025=@(c,interpolated_values) -0.75*exp(c(1)*interpolated_values(:,2).^2.*interpolated_values(:,4)+c(2)*interpolated_values(:,3).^2.*interpolated_values(:,4)) +1 ;
%%
% estimating G(alpha,beta) for N=25 
[c1_n25,resnorm1_n25] = lsqcurvefit(fun1, [-1 -1] ,interpolated_values_n25,interpolated_values_n25(:,5)); 
[c1_n50,resnorm1_n50] = lsqcurvefit(fun1, [-1 -1] ,interpolated_values_n50,interpolated_values_n50(:,5));
[c1_n100,resnorm1_n100] = lsqcurvefit(fun1, [-1 -1] ,interpolated_values_n100,interpolated_values_n100(:,5)); 
[c1_n150,resnorm1_n150] = lsqcurvefit(fun1, [-1 -1] ,interpolated_values_n150,interpolated_values_n150(:,5)); 
[c1_n250,resnorm1_n250] = lsqcurvefit(fun1, [-1 -1] ,interpolated_values_n250,interpolated_values_n250(:,5)); 
[c1_n500,resnorm1_n500] = lsqcurvefit(fun1, [-1 -1] ,interpolated_values_n500,interpolated_values_n500(:,5)); 
c1_n=[25 c1_n25;50 c1_n50; 100 c1_n100;150 c1_n150 ;250 c1_n250;500 c1_n500];
%%
[c2_n ,resnorm2_n] = lsqcurvefit(fun2, [-1 -1] ,interpolated_values,interpolated_values(:,5)); 
%%
[c2_n_lowerbound_025 ,resnorm2_n_lowerbound_025] = lsqcurvefit(fun2_lowerbound_025, [-1 -1] ,interpolated_values,interpolated_values(:,5)); 
