% double shift qr
clc;
clear;

pkg load signal

n = 8;
D = dctmtx(n); % we use DCT matrix since it has complex eigenvalues

%% double shift qr
function y = doubleshiftqr(A, max_iter = 20)
  mu = zeros(2,1);
  n = size(A,1);
  ei = zeros(n,1);
  if n == 1
    ei = A(1,1);
  else
    I = diag(ones(n,1));
    for i = 1:max_iter
      % solving characteristic equation to get mu1, mu2
      mu = roots([1 -(A(end-1,end-1)+A(end,end)) (A(end-1,end-1)*A(end,end)-A(end-1,end)*A(end,end-1))]);
      [Q, R] = qr((A-mu(1)*I)*(A-mu(2)*I));
      A = Q'*A*Q;
    endfor
  endif
  
  i = 1;
  while i < n
    if abs(A(i+1,i)) < 1e-16
      ei(i) = A(i,i);
      i = i + 1;
    else
      temp = zeros(2,2); lambda = zeros(2,1);
      temp(1,1) = A(i,i); temp(1,2) = A(i,i+1); temp(2,1) = A(i+1,i); temp(2,2) = A(i+1,i+1);
      lambda = roots([1 -trace(temp) det(temp)]);
      ei(i) = lambda(1);
      ei(i+1) = lambda(2);
      i = i + 2;
    endif
  endwhile
  if abs(A(n,n-1)) < 1e-16
    ei(n) = A(n,n);
  endif
  A
  y = sort(ei);
endfunction

res = doubleshiftqr(D, 750)
error = norm((sort(real(res)) - sort(real(eig(D)))) + (sort(imag(res)) - sort(imag(eig(D)))))