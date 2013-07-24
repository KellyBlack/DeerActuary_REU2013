fp=fopen('Xem.csv','w');
fprintf(fp,'Xem\n');

for j=1:500
    [X, Xem]=TestX(1, 100, .05);
    fprintf(fp,'%f\n',Xem(101));
end

fclose(fp);