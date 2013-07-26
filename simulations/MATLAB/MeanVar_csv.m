function [  ] = MeanVar_csv( alpha, alpha0, dalpha, gamma, gamma0, dgamma, P, P0, dP )
%Wraps MeanVar_DeerInsSim and records means and variances for different
%values of alhpa, gamma, and P (our unknown parameters)
%   Detailed explanation goes here

%{
do this at least once before the function to create header
fp = fopen('meanvarEM.csv','w');
fprintf(fp,'alpha, gamma, P, Xmean, Xvar, Mmean, Mvar\n');
fclose(fp);
%}

tic

for i=alpha0:dalpha:alpha
    for j=gamma0:dgamma:gamma
        for k=P0:dP:P
            fp=fopen('meanvarEM.csv','a');
            [Xmean, Xvar, Mmean, Mvar] = MeanVar_DeerInsSim(i, j, k);
            fprintf(fp,'%f, %f, %f, %f, %f, %f, %f\n',i, j, k, Xmean, Xvar, Mmean, Mvar);
            fclose(fp);
        end
    end
end

toc