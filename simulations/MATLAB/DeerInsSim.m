function [ Xtrue, Mem ] = DeerInsSim(  T, N, r1, h, F, alpha, rho, beta, P, gamma, g )
%Calls X and Mem with same W
%   Detailed explanation goes here

dt = T/N;
%s = (0:dt:T); %time
dW = sqrt(dt)*randn(1, N+1);       
W = [0, cumsum(dW(1:N))]; % W is initially 0

[ Xtrue ] = DeerPop_Beth(T, N, r1, h, F, alpha, dW, W);
[ Mem ] = InsPayout(T, N, r1, h, F, alpha,  rho, beta, P, gamma, g, dW, W);



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
