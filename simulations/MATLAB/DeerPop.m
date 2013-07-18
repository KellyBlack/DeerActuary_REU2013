function [ Z, X ] = DeerPop( T, N, r1, h, F, alpha )
% Solution to our deer population Logistic SDE with harvest
%   Detailed explanation goes here

% rescale the parameters for the new harvest rate
r2=r1-h;              % r tilde
f=(r2/r1)*F;          % f tilde
a=(alpha^2)-r2;       % a from linearization

b=(.5*(alpha^2))-a;   % quantity used a lot...

g0= 1-(r2/b); % constant so that Z is initially 1
 
dt=sqrt(T/N);
dW = dt*randn(N+1,1);
%dW(1) = 0.0;         % W is initially 0
W = [0; cumsum(dW(1:N))]; 

ds=(T/N);
s=(0:ds:T)';    % time 

Q = cumsum(exp(b*s+alpha*W).*dW);   %integral part of analytic solution
Z = (r2/b) + g0*exp((-b*s)-(alpha*W)) + exp((-b*s)-(alpha*W)).*((-alpha*r2)/b).*Q; %analytic solution of Z
X =(f./Z); % transform back to X




