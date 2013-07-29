function [  ] = MeanVarMil_csv( alpha, alpha0, dalpha, gamma, gamma0, dgamma, P, P0, dP, iterations )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

%{
do this at least once before the function to create header
fp = fopen('meanvarmil100.csv','w');
fprintf(fp,'alpha, gamma, P, Xmean, Xvar, Mmean, Mvar\n');
fclose(fp);
%}

tic

for i=alpha0:dalpha:alpha
    for j=gamma0:dgamma:gamma
        for k=P0:dP:P
            fp=fopen('meanvarmil_25000_500.csv','a');
            [Xmean, Xvar, Mmean, Mvar] = MeanVar_DeerInsSim_Mil(i, j, k, iterations);
            fprintf(fp,'%f, %f, %f, %f, %f, %f, %f\n',i, j, k, Xmean, Xvar, Mmean, Mvar);
            fclose(fp);
        end
    end
end

toc