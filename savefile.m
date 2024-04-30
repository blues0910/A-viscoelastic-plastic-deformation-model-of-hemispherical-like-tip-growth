function savefile(s,d)
fid = fopen(s,"w");
[n,m] = size(d);
if n < m
    d = d';
end
if min(n,m)==2
    for i = 1:max(n,m)
        fprintf(fid,'%f %f',d(i,1),d(i,2));
        fprintf(fid,'\n');
    end
elseif min(n,m)==1
    for i = 1:max(n,m)
        fprintf(fid,'%f',d(i));
        fprintf(fid,'\n');
    end
end
fclose(fid);
end