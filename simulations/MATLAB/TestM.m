function [ M, Mem ] = TestM( T, N, alpha, rho )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

dt = T/N;
s = (0:dt:T);
dW = sqrt(dt)*randn(1, N+1);
W = [0, cumsum(dW(1:N))];

[ X, Xem ] = TestX( T, N, alpha);

m0=0;
Mem = zeros(1, N+1);
Mem(1) = m0;
Mtemp = m0;

for i=1:N
  Winc = dW(i);
  Mtemp = Mtemp + dt * ((rho * Mtemp) - ((rho/alpha) * X(i))) + Winc;
  Mem(i + 1) = Mtemp; 
end

M = m0*exp(rho*s) + W;