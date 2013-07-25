function [ Xtrue, Mmil ] = DeerInsSimMil(  T, N, r1, h, F, alpha, rho, beta, P, gamma, g )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

dt = T/N;
s = (0:dt:T);

[ Xtrue, Xmil ] = DeerPopMil(T, N, r1, h, F, alpha);
[ Mmil ] = InsPayoutMil(T, N, r1, h, F, alpha,  rho, beta, P, gamma, g);



%{

clf;

%Deer Population
subplot(1, 2, 1);
plot(s, Xtrue, 'g-', s, Xmil, 'rx')
title('Deer Population');
legend('Xtrue', 'Xmil')

%Insurance Payout
subplot(1, 2, 2);
plot(s, Mmil,'rx')
title('Insurance Payout');
legend('Mem')

%}
