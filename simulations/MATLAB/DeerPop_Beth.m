function [ X ] = DeerPop_Beth( T, N, r1, h, F, alpha, dW, W )
% Solution to our deer population Logistic SDE with harvest and EM
% approximation
%   Detailed explanation goes here

%Beth's code
% rescale the parameters for the new harvest rate
rtilde=r1-h;              % r tilde
ftilde=(rtilde/r1)*F;          % f tilde
a=(alpha^2)-rtilde;       % a from linearization

b=(.5*(alpha^2))-a;   % quantity used a lot...

g0= 1-(rtilde/b); % constant so that Z is initially 1
 
dt=T/N;
%dW = sqrt(dt)*randn(1, N+1);       
%W = [0, cumsum(dW(1:N))]; % W is initially 0
s=(0:dt:T);    % time 

Q = cumsum(exp(b*s+alpha*W).*dW);   %integral part of analytic solution
Z = (rtilde/b) + g0*exp((-b*s)-(alpha*W)) + exp((-b*s)-(alpha*W)).*((-alpha*rtilde)/b).*Q; %analytic solution of Z
X = (ftilde./Z); % transform back to X

%clf;
%plot(s, X, 'g-')
%hold on


%Candace's code for EM approx of X
% x0 = ftilde;
% 
% Xem = zeros(1,N+1); %initialize EM vector
% Xem(1) = x0;    %setting the initial condition
% Xtemp_em = x0;
% 
% 
% %EM Approximation
% for j = 1:N
%    Winc = dW(j);
%    Xtemp_em = Xtemp_em + ds * rtilde * Xtemp_em *(1-Xtemp_em/ftilde) + alpha * Xtemp_em * Winc;
%    Xem(j+1) = Xtemp_em;
% end

%plot (s, Xem, 'r*')    
%legend('Xtrue','Xem') 