% Power method
clc;
clear;

function [vec,value]=power(start,A,tolerr)
  dd=1;
  x=start;
  n_=0;
  n=10;
  i = 0;
  while dd>tolerr
    n_ = n;
    y=A*x;
    dd=abs(norm(x)-n);
    n=norm(x);
    i = i+1
    x=y/n;
    norm(n-n_)/norm(n)
  end
  vec=x;
  value=n;
endfunction
A1 = [11 7 -4
7 11 4
-4 4 10];
A2 = [2 2
2 -1];
A3 = [47 32 8
32 -1 -16
8 -16 59];
R = rosser;
[x y] = power([1;1;1;1;1;1;1;1],R,1e-9)
eigvec = (A3-y*eye(3))(:,3)
