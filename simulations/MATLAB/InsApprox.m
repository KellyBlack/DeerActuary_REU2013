function [ Mem, Mmil ] = InsApprox( T, N, rho, beta, P, gamma )
%UNTITLED5 Summary of this function goes here
%   T time
% N Steps
% rho insurance rate thing (capital R)
% beta claims rate
% P payout
% gamma noise coefficient

%finer partition 
R = 4;
L = R*N;

dt = T/N;
Dt = dt/R;
dW = sqrt(dt)*randn(L+1,1);
dW(1)=0.0;

m0 = 225000000;

Mem = zeros(L+1, 1); %initialize EM vector
Mem(1) = m0;    %setting the initial condition
Mtemp_em = m0;

Mmil = zeros(L+1, 1); %initilize Mil vector
Mmil(1) = m0;    %setting the initial condition
Mtemp_mil = m0;

[X, Z] = DeerPop(10, L, 30, 20, 3000000, 1); %call X

%EM and Milstein Approximations
for j = 1:L
   Winc = dW(j);
   Mtemp_em = Mtemp_em + Dt * ((rho * Mtemp_em) -(beta * X(j)) + P) - (gamma * X(j)) * Winc;
   Mem(j+1) = Mtemp_em;
   Mtemp_mil = Mtemp_mil + Dt * ((rho * Mtemp_mil) -(beta * X(j)) + P) - (gamma * X(j)) * Winc + 0.5*gamma * gamma * X(j) *(Winc*Winc-Dt);
   Mmil(j+1) = Mtemp_mil;
end

time = 0:Dt:T;     %time for EM and Mil Approximations
plot (time, Mem, 'r*',time, Mmil, 'b*')    
legend('Mem','Mmil')
