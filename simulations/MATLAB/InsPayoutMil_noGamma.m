function [ Mmil ] = InsPayoutMil_noGamma( T, N, r1, h, F, alpha, rho, beta, P, g, dW )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

rtilde=r1-h;
ftilde=(rtilde/r1)*F;

 dt=T/N;
% dW = dt*randn(1, N+1);

m0 = (beta*ftilde - P)/(rho-g);
Mmil = zeros(1, N+1);
Mmil(1) = m0;
Mtemp = m0;

[ Xtrue ] = DeerPopMil( T, N, r1, h, F, alpha, dW );

for i=1:N
    Winc = dW(i);
    Mtemp = Mtemp + dt*(rho*Mtemp - beta*Xtrue(i)+P) - Winc*beta*Xtrue(i) - .5*beta*alpha*Xtrue(i)*(Winc*Winc-dt);
    Mmil(i+1) = Mtemp;
end

% s=(0:dt:T);
% plot(s, Mmil, 'bx')
% legend('Mmil')