% QR method

clc;
clear;

pkg load signal

n = 8;
D = dctmtx(n); % we use DCT matrix since it has complex eigenvalues

max_iter = 50000;
tol_err = 1e-12;

C = D;
n = size(C,1);
ei = zeros(n,1);
I = eye(n);

for i = 1:max_iter
  [Q, R] = qr(C);
  C = R*Q;
endfor

i = 1;
while i < n
  if abs(C(i+1,i)) < 1e-15
    ei(i) = C(i,i);
    i = i + 1;
  else
    temp = zeros(2,2); lambda = zeros(2,1);
    temp(1,1) = C(i,i); temp(1,2) = C(i,i+1); temp(2,1) = C(i+1,i); temp(2,2) = C(i+1,i+1);
    lambda = roots([1 -trace(temp) det(temp)]);
    ei(i) = lambda(1);
    ei(i+1) = lambda(2);
    i = i + 2;
  endif
endwhile
if C(n,n-1) == 0
  ei(n) = C(n,n);
endif

res = ei
error = norm((sort(real(res)) - sort(real(eig(D)))) + (sort(imag(res)) - sort(imag(eig(D)))))