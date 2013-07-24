fp=fopen('Mem.csv','w');
fprintf(fp,'Mem\n');

for j=1:1000
    [M, Mem]=TestM(1, 100, .05, .1);
    fprintf(fp,'%f\n',Mem(101));
end

fclose(fp);