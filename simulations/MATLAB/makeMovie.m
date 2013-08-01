
T = 2.0;
N = 10000; %start with 10000 then use 25000
r1 = 1.7;
h = .16;
F = 28000;
rho = .04; %log(1+rate)
beta = 9; %.003*3000

% Set up the parameters for the deer simulation
alpha = 0.15;
gamma = 1.20;
P     = 277500.0;

rtilde=r1-h;
ftilde=(rtilde/r1)*F;
a=(alpha^2)-rtilde;      
b=(.5*(alpha^2))-a;
g0= 1-(rtilde/b);

% Variables for the bond fund

g = .05;
m0 = (beta*ftilde - P)/(rho-g);
rho = .004; %log(1+rate)
beta = 9; %.003*3000
m = m0;

dt=T/N;
sdt = sqrt(dt);

clf;
h = figure(1);

%Deer Population on 1
subplot(2, 1, 1);
axis([0 T ftilde*0.7 ftilde*1.3]);
plot([0 T],[ftilde ftilde],'k--');
title('A Simulation Of The Deer Population and Resulting Fund Balance');
ylabel('Deer Population');
xlabel('time');
hold on;

%Bond funds on 2
subplot(2, 1, 2);
axis([0 T 0.98*m0 1.04*m0]);
ylabel('Bond Fund Balance');
xlabel('time');
hold on;

% Set up the movie buffers
skipFrame = 8;
%set(gca,'NextPlot','replaceChildren');
frame(round(N/skipFrame/2)) = struct('cdata',[],'colormap',[]);
frameNumber = 1;


fprintf('Beginning simulation\n');
z = 1.0;
X = ftilde;
W = 0.0;
Q = 0.0;
for j=0:N,
    dW = sdt*randn(1,1);
    s = j*dt;
    
    % Update the deer population
    Q = Q + exp(b*s+alpha*W).*dW;   %integral part of analytic solution
    Z = (rtilde/b) + g0*exp((-b*s)-(alpha*W)) + ...
        exp((-b*s)-(alpha*W)).*((-alpha*rtilde)/b).*Q; %analytic solution of Z
    prevX = X;
    X = (ftilde./Z); % transform back to X
    
    % Update the bond fund
    m = m + dt*(rho*m - beta*prevX+P) - dW*gamma*prevX - .5*gamma*alpha*prevX*(dW*dW-dt);
    W = W + dW;

    if(mod(j,skipFrame)==0)
        % Plot the updates
        subplot(2, 1, 1);
        plot(s,X,'rd','MarkerSize',2);
        axis([0 T ftilde*0.7 ftilde*1.3]);
        subplot(2, 1, 2);
        plot(s,m,'bo','MarkerSize',2);
        drawnow();
        if(mod(j,2*skipFrame)==0)
            frame(frameNumber) = getframe(h);
            frameNumber = frameNumber + 1;
        end
    end
    
end
