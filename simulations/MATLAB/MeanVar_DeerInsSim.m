function [ Xmean, Xvar, Mmean, Mvar ] = MeanVar_DeerInsSim( alpha, gamma, P )
%Wraps DeerInsSim to find mean and variance of X and Mem
%   Detailed explanation goes here

iterations = 500;

T = 10;
N = 25000; %start with 10000 then use 25000
r1 = 1.7;
h = .16;
F = 28000;
rho = .004; %log(1+rate)
beta = 9; %.003*3000
g = .005;

Xvec = zeros(1, iterations); %initialize X vector for final positions
Mvec = zeros(1, iterations); %initialize M vector for final positions

for i=1:iterations
    [Xtrue, Mem] = DeerInsSim( T, N, r1, h, F, alpha, rho, beta, P, gamma, g );
    Xvec(1,i) = Xtrue(N+1);
    Mvec(1,i) = Mem(N+1);
end

Xmean = mean(Xvec);
Xvar = var(Xvec);
Mmean = mean(Mvec);
Mvar = var(Mvec);


