function [ X, Xem, M, Mem ] = Testing( T, N, alpha, rho )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
dt = T/N;
s = (0:dt:T);

clf;

%X
subplot(1, 2, 1);
[X, Xem] = TestX( T, N, alpha );
plot(s, X, 'g-')
hold on
plot(s, Xem, 'r*')
title('X');
legend('X', 'Xem')

%M
subplot(1, 2, 2);
[M, Mem] = TestM(T, N, alpha, rho);
plot(s, M,'g-')
hold on
plot( s, Mem, 'r*')
title('M');
legend( 'M', 'Mem')
