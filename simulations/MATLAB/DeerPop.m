function [ Ztrue, Xtrue ] = DeerPop( T, N, r1, h, F, alpha )
% Solution to our deer population Logistic SDE with harvest
%   Detailed explanation goes here
%randn('state',100)

r2=r1-h;
f=(r2/r1)*F;
a=(alpha^2)-r2;  
b=(.5*(alpha^2))-a;

g0= 1-(r2/b); 
% problem parameters
 
dt=sqrt(T/N);
dW = dt*randn(1,N);        
W = cumsum(dW); 

ds=T/N;
s=ds:ds:T;    %time

for i=1:N
    Q(i)=cumsum(exp(b*s(i) + alpha*W(i))*dW(i));    
end


Ztrue = (r2/b) + g0*exp((-b*T)-(alpha*W)) + exp((-b*T)-(alpha*W))*((-alpha*r2)/b)*Q(i);

Xtrue=f + (f./Ztrue);


end

