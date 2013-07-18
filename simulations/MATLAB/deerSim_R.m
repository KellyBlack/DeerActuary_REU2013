function [ Ztrue, Xtrue, Xem, Xmil  ] = deerSim_R( T, N, r1, h, F, alpha )
%Simulation of Deer Population
%   Uses Xtrue (better approximation), Euler-Maruyama, Milstein 

%define rtilde and ftilde with input
rtilde = r1-h;
ftilde = ((r1-h)/r1)*F;

%finer partition for true data
R = 4;
L = R*N;

dt = T/N;
Dt = dt/R;
dW = sqrt(dt)*randn(1,L+1);
dW(1)=0.0;
W = cumsum(dW);

time1 = 0:dt:T;     %time for z

x0 = ftilde;
z0 = ftilde/x0;
a = alpha^2 -rtilde;
b = alpha;
g0 = 1-rtilde/(.5*(b*b)-a);
phi = exp((a - .5*b*b) * time1 - b*W(1:R:(L+1))); %exp in integral

integral = zeros(1,N+1); %initializes the integral vector
integral(1) = 0;
s = 0;
ds = dt;

%Evaluates the integral in z
for i=1:N
    integral(i+1) = integral(i) + (1.0/phi(i)) * sum(dW(R*(i-1)+1:R*i));
    s = i*ds;
end

%evaluates Z in its entirety
Ztrue = rtilde/(.5*b*b - a) + g0 * phi - b * rtilde/(.5*b*b - a) * phi .* integral;

%with Ztrue, Xtrue can be calculated
Xtrue = ftilde./Ztrue;

clf;
plot(time1, Xtrue, 'g-')
hold on

Xem = zeros(1,L+1); %initialize EM vector
Xem(1) = x0;    %setting the initial condition
Xtemp_em = x0;

Xmil = zeros(1,L+1); %initilize Mil vector
Xmil(1) = x0;    %setting the initial condition
Xtemp_mil = x0;

%EM and Milstein Approximations
for j = 1:L
   Winc = dW(j);
   Xtemp_em = Xtemp_em + Dt * rtilde * Xtemp_em *(1-Xtemp_em/ftilde) + alpha * Xtemp_em * Winc;
   Xem(j+1) = Xtemp_em;
   Xtemp_mil = Xtemp_mil + Dt * rtilde * Xtemp_mil * (1-Xtemp_mil/ftilde) + alpha *Xtemp_mil* Winc + 0.5*alpha * alpha * Xtemp_mil *(Winc*Winc-Dt);
   Xmil(j+1) = Xtemp_mil;
end

time2 = 0:Dt:T;     %time for EM and Mil
plot (time2, Xem, 'r*',time2, Xmil, 'b*')    
legend('Xtrue','Xem','Xmil')


