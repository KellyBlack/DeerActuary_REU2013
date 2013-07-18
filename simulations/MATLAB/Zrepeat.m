fp=fopen('Z.csv','w');
fprintf(fp,'Z\n');

for i=1:10000
    [Z, X] = DeerPop(10, 1000, 30, 20, 3000000, 1);
     fprintf(fp,'%f\n',Z(1000+1)); %records last position Z
end

fclose(fp); 