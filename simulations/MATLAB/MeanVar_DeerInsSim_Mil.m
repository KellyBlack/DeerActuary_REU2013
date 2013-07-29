function [ Xmean, Xvar, Mmean, Mvar ] = MeanVar_DeerInsSim_Mil( alpha, gamma, P, iterations )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here


T = 10;
N = 25000; %start with 10000 then use 25000
r1 = 1.7;
h = .16;
F = 28000;
rho = .004; %log(1+rate)
beta = 9; %.003*3000
g = .005;

dt = sqrt(T/N);

Xvec = zeros(1, iterations);
Mvec = zeros(1, iterations);

for i=1:iterations
    dW = dt*randn(1, N+1);
    [Xtrue, Mmil] = DeerInsSimMil( T, N, r1, h, F, alpha, rho, beta, P, gamma, g, dW );
    Xvec(1,i) = Xtrue(N+1);
    Mvec(1,i) = Mmil(N+1);
end

Xmean = mean(Xvec);
Xvar = var(Xvec);
Mmean = mean(Mvec);
Mvar = var(Mvec);


