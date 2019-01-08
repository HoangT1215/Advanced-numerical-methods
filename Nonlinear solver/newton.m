% 2-1 preclass: bisection method

function [ sol ] = NewtonMethod( f,f0,tol )
syms x;
syms f(x) = matlabFunction(f);
f_dash(x)=diff(f(x));
sol=f0;
ratio=double(((f(f0))/(f_dash(f0))));
while double(abs(ratio))>= tol;
    sol=double(sol-ratio);
    ratio=double((f(sol))/(f_dash(sol)));
end
end

NewtonMethod(x^2 - 2, 1, 1e-6)