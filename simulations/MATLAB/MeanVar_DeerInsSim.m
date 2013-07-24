function [ Xmean, Xvar, Mmean, Mvar ] = MeanVar_DeerInsSim( alpha, gamma, P )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

%Problem Parameters
T = 10;
N = 10000; %start with 10000 then use 25000
r1 = 1.7;
h = .16;
F = 27356;
rho = .004; %log(1+rate)
beta = 200; %
g = .04;

for i=1:10
    [Xtrue, Mem] = DeerInsSim( T, N, r1, h, F, alpha, rho, beta, P, gamma, g );
end

Xmean = mean(Xtrue);
Xvar = var(Xtrue);
Mmean = mean(Mem);
Mvar = var(Mem);


