%Linton's Insurance Payout Code from R
T =10;
N = 3600;
r1 = 30; 
h = 20;
F = 3000000;
alpha = 1;

rtilde = r1-h;
ftilde = (rtilde/r1)*F;

[ Ztrue , Xtrue ] = DeerPop_Beth( T, N, r1, h, F, alpha );

L = N;
dt = T/N;
s = (0:dt:T);

dW = sqrt(dt)*randn(1,L+1);

I = 0.0005;
beta = 0.1;
P = 10000;
gamma = 0.1;
g = 0.0001;

m0 = (beta * ftilde - P)/(I - g);
Mem = zeros(1, N+1);
Mmil = zeros(1, N+1);
Mem(1) = m0;
Mmil(1) = m0;
Mtemp1 = m0;
Mtemp2 = m0;

for i=1:L
  %Winc = sum(dw[(R * (i - 1) + 1):(R * i)]);
  Winc = dW(i);
  Mtemp1 = Mtemp1 + dt * (I * Mtemp1 - beta * Xtrue(i) + P) - Winc * gamma * Xtrue(i);
  Mem(i + 1) = Mtemp1;
  Mtemp2 = Mtemp2 + dt * (I * Mtemp2 - beta * Xtrue(i) + P) - Winc * gamma * Xtrue(i) + 0.5 * gamma * gamma * Xtrue(i) * (Winc * Winc - dt);
  Mmil(i + 1) = Mtemp2;
end


clf;
plot(s, Mem, 'r*', s, Mmil, 'b*')




