function [  ] = MeanVar_DeerInsSim_Mil_noGamma( alpha, P, iterations )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
% 
% fp = fopen('MeanVar_noGamma.csv','w');
% fprintf(fp,'X, M, Z\n');
% fclose(fp);

T = 10;
N = 25000; %start with 10000 then use 25000
r1 = 1.7;
h = .16;
F = 28000;
rho = .004; %log(1+rate)
beta = 17; %.0057*3000
g = .005;

dt = T/N;
Dt=sqrt(dt);

% Xvec = zeros(1, iterations);
% Mvec = zeros(1, iterations);
tic

for i=1:iterations
    dW = Dt*randn(1, N+1);
    
        fp=fopen('MeanVar_noGamma.csv','a');
        
    [Xtrue, Mmil] = DeerInsSimMil_noGamma( T, N, r1, h, F, alpha, rho, beta, P, g, dW );
    
        fprintf(fp,'%f, %f\n', Xtrue(N+1), Mmil(N+1));
        fclose(fp);
%     Xvec(1,i) = Xtrue(N+1);
%     Mvec(1,i) = Mmil(N+1);
end
toc

% Xmean = mean(Xvec);
% Xvar = var(Xvec);
% Mmean = mean(Mvec);
% Mvar = var(Mvec);


