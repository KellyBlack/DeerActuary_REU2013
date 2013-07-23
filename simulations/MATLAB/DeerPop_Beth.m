function [X , Xem] = DeerPop_Beth( T, N, r1, h, F, alpha )
% Solution to our deer population Logistic SDE with harvest
%   Detailed explanation goes here

%Beth's code
% rescale the parameters for the new harvest rate
rtilde=r1-h;              % r tilde
ftilde=(rtilde/r1)*F;          % f tilde
a=(alpha^2)-rtilde;       % a from linearization

b=(.5*(alpha^2))-a;   % quantity used a lot...

g0= 1-(rtilde/b); % constant so that Z is initially 1
 
dt=sqrt(T/N);
dW = dt*randn(1, N+1);
%dW(1) = 0.0;         % W is initially 0
W = [0, cumsum(dW(1:N))]; 

ds=T/N;
s=(0:ds:T);    % time 

Q = cumsum(exp(b*s+alpha*W).*dW);   %integral part of analytic solution
Z = (rtilde/b) + g0*exp((-b*s)-(alpha*W)) + exp((-b*s)-(alpha*W)).*((-alpha*rtilde)/b).*Q; %analytic solution of Z
X = (ftilde./Z); % transform back to X

%clf;
%plot(s, X, 'g-')
%hold on
%%{

%Candace's code
x0 = ftilde;

Xem = zeros(1,N+1); %initialize EM vector
Xem(1) = x0;    %setting the initial condition
Xtemp_em = x0;

%Xmil = zeros(1,N+1); %initialize Mil vector
%Xmil(1) = x0;    %setting the initial condition
%Xtemp_mil = x0;

%EM and Milstein Approximations
for j = 1:N
   Winc = dW(j);
   Xtemp_em = Xtemp_em + ds * rtilde * Xtemp_em *(1-Xtemp_em/ftilde) + alpha * Xtemp_em * Winc;
   Xem(j+1) = Xtemp_em;
   %Xtemp_mil = Xtemp_mil + ds * rtilde * Xtemp_mil * (1-Xtemp_mil/ftilde) + alpha *Xtemp_mil* Winc + 0.5*alpha * alpha * Xtemp_mil *(Winc*Winc-ds);
   %Xmil(j+1) = Xtemp_mil;
end

%plot (s, Xem, 'r*',s, Xmil, 'b*')    
%legend('Xtrue','Xem','Xmil')


%}