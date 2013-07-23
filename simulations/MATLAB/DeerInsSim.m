function [ Xtrue, Mem ] = DeerInsSim(  T, N, r1, h, F, alpha, rho, beta, P, gamma, g )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

dt = T/N;
s = (0:dt:T);

[ Xtrue, Xem ] = DeerPop_Beth(T, N, r1, h, F, alpha);
[ Mem ] = InsPayout(T, N, r1, h, F, alpha,  rho, beta, P, gamma, g);



%{

clf;

%Deer Population
subplot(1, 2, 1);
plot(s, Xtrue, 'g-', s, Xem, 'rx')
%hold on
%plot(s, Xmil, 'b*')
title('Deer Population');
legend('Xtrue', 'Xem')

%Insurance Payout
subplot(1, 2, 2);
plot(s, Mem,'rx')
%hold on
%plot( s, Mmil, 'b*')
title('Insurance Payout');
legend( 'Mem')

%}
