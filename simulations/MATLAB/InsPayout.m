function [ Mem ] = InsPayout( T, N, r1, h, F, alpha, rho, beta, P, gamma, g )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%Linton's Insurance Payout Code from R

rtilde = r1-h;
ftilde = (rtilde/r1)*F;

[ Xtrue, Xem ] = DeerPop_Beth( T, N, r1, h, F, alpha );

L = N;
dt = T/N;
s = (0:dt:T);

dW = sqrt(dt)*randn(1,L+1);

%rho = 0.0005; %bond rate
%beta = 0.1; %claim rate
%P = 10000; %principal premium
%gamma = 0.5; %noise parameter
%g = 0.0001; %rate of profit for the company

m0 = (beta * ftilde - P)/(rho - g);
Mem = zeros(1, N+1);
%Mmil = zeros(1, N+1);
Mem(1) = m0;
%Mmil(1) = m0;
Mtemp1 = m0;
%Mtemp2 = m0;

for i=1:L
  %Winc = sum(dw[(R * (i - 1) + 1):(R * i)]);
  Winc = dW(i);
  Mtemp1 = Mtemp1 + dt * (rho * Mtemp1 - beta * Xtrue(i) + P) - Winc * gamma * Xtrue(i);
  Mem(i + 1) = Mtemp1;
  %Mtemp2 = Mtemp2 + dt * (rho * Mtemp2 - beta * Xtrue(i) + P) - Winc * gamma * Xtrue(i) + 0.5 * gamma * gamma * Xtrue(i) * (Winc * Winc - dt);
  %Mmil(i + 1) = Mtemp2;
end


%clf;
%plot(s, Mem, 'r*', s, Mmil, 'b*')
%legend('Mem', 'Mmil');


