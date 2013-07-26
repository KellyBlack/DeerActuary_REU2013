function [ Mem ] = InsPayout( T, N, r1, h, F, alpha, rho, beta, P, gamma, g, dW, W )
%EM Approximation for M
%   Detailed explanation goes here


rtilde = r1-h;
ftilde = (rtilde/r1)*F;

[ Xtrue ] = DeerPop_Beth( T, N, r1, h, F, alpha, dW, W );

dt = T/N;
%s = (0:dt:T);
%dW = sqrt(dt)*randn(1,N+1);

m0 = (beta * ftilde - P)/(rho - g); %initial M
Mem = zeros(1, N+1); %create Mem vector
Mem(1) = m0; %initialize vector
Mtemp1 = m0;


for i=1:N
  Winc = dW(i);
  Mtemp1 = Mtemp1 + dt * (rho * Mtemp1 - beta * Xtrue(i) + P) - Winc * gamma * Xtrue(i);
  Mem(i + 1) = Mtemp1;
end


%clf;
%plot(s, Mem, 'r*')
%legend('Mem');


