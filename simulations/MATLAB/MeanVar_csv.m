function [  ] = MeanVar_csv( alpha, alpha0, dalpha, gamma, gamma0, dgamma, P, P0, dP )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

%{
do this at least once before the function to create header
fp = fopen('name.csv','w');
fprintf(fp,'alpha, gamma, P, Xmean, Xvar, Mmean, Mvar\n')
fclose(fp);
%}


for i=alpha0:dalpha:alpha
    for j=gamma0:dgamma:gamma
        for k=P0:dP:P
            fp=fopen('meanvar.csv','a');
            [Xmean, Xvar, Mmean, Mvar] = MeanVar_DeerInsSim(i, j, k);
            fprintf(fp,'%f, %f, %f, %f, %f, %f, %f\n',alpha, gamma, P, Xmean, Xvar, Mmean, Mvar);
            fclose(fp);
        end
    end
end