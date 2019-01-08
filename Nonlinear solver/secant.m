% Secant method
clc;
clear;

function y = f(x)
  y = x^2-1-exp(x);
endfunction

x0 = 2; % initial value
x = 3;
x_ = x0;
tolerr = 1e-6;
relerr = 0;
c = (1+sqrt(5))/2
n = 0;

while abs(f(x)) > tolerr
  temp = x;
  x = x - f(x)*(x-x_)/(f(x)-f(x_));
  x_ = temp;
  relerr = abs(x-x_)/x
  n = n + 1;
endwhile
