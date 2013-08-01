function [ Xtrue, Mmil ] = DeerInsSimMil_noGamma(  T, N, r1, h, F, alpha, rho, beta, P, g, dW )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% dt = T/N;
% s = (0:dt:T);

[ Xtrue ] = DeerPopMil(T, N, r1, h, F, alpha, dW);
[ Mmil ] = InsPayoutMil_noGamma(T, N, r1, h, F, alpha, rho, beta, P, g, dW);





% clf;
% 
% %Deer Population
% subplot(2, 1, 1);
% plot(s, Xtrue, 'g-')
% title('Deer Population');
% legend('Xtrue')
% xlabel('t')
% 
% %Insurance Payout
% subplot(2, 1, 2);
% plot(s, Mmil,'rx')
% title('Insurance Payout');
% legend('Mmil')
% xlabel('t')


