function [ X, Xem ] = TestX( T, N, alpha)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

dt=sqrt(T/N);
dW = dt*randn(1, N+1);
W = [0, cumsum(dW(1:N))]; 

%ds=T/N;
%s=(0:ds:T);    % time 

x0=0;
Xem = zeros(1,N+1); %initialize EM vector
Xem(1) = x0;    %setting the initial condition
Xtemp_em = x0;



%EM Approximation
for j = 1:N
   Winc = dW(j);
   Xtemp_em = Xtemp_em + alpha * Winc;
   Xem(j+1) = Xtemp_em;
end

X=alpha.*W;