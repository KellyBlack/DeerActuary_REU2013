function [ ] = Simple( T, N, M )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

fp=fopen('Simple.csv','w');
fprintf(fp,'Z\n');

dt=T/N;
for i=1:M
    %time=1:dt:T;
    dW = sqrt(dt)*randn(N+1,1);
    W = [0; cumsum(dW(1:N))]; 
    Z=cumsum(W.*dW);
    
    %dW(1) = 0.0;         % W is initially 0
    %W1 = cumsum(dW);
    %Z1=cumsum(W1.*dW);
    
    fprintf(fp,'%f\n',Z(N+1)); %records last position Z
    %fprintf(fp,'%f,%f\n',Z(N+1),Z1(N+1)); %records last position Z
end

fclose(fp);
