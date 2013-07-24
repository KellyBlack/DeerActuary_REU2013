function [ Mem, Mmil ] = InsApprox( T, N, rho, beta, P, g, gamma )
%UNTITLED5 Summary of this function goes here
%   T time
% N Steps
% rho bond rate (capital R)
% beta claim rate 
% P principle premium
% gamma noise coefficient

%finer partition 
%R = 4;
%L = R*N;

dt = T/N;
%Dt = dt/R;
dW = sqrt(dt)*randn(N+1,1);


f=1000000;
m0 =((beta * f) - P)/(rho - g); %g rate or profit

Mem = zeros(N+1, 1); %initialize EM vector
Mem(1) = m0;    %setting the initial condition
Mtemp_em = m0;

Mmil = zeros(N+1, 1); %initilize Mil vector
Mmil(1) = m0;    %setting the initial condition
Mtemp_mil = m0;

[X, Z] = DeerPop(T, N, 30, 20, 3000000, 1); %call X

%EM and Milstein Approximations
for j = 1:N
   Winc = dW(j);
   Mtemp_em = Mtemp_em + (dt * ((rho * Mtemp_em) -(beta * X(j)) + P)) - ((gamma * X(j)) * Winc);
   Mem(j+1) = Mtemp_em;
   Mtemp_mil = Mtemp_mil + (dt * ((rho * Mtemp_mil) -(beta * X(j)) + P)) - ((gamma * X(j)) * Winc) + (0.5*gamma * gamma * X(j) *(Winc*Winc-dt));
   Mmil(j+1) = Mtemp_mil;
end

time = 0:dt:T;     %time for EM and Mil Approximations
plot (time, Mem, 'r*',time, Mmil, 'bo')    
legend('Mem','Mmil')
