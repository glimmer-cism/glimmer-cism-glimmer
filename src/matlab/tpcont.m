function [data] = tpcont(file,year,varb);

id = fopen(file,'r','ieee-be');
a = fread(id,1,'int');

if (a > 100)
  fclose(id);
  id = fopen(file,'r','ieee-le');
else
  frewind(id);
end

err = 0;
[t,s,x,y] = tphead(id);

while (year ~= t) | (strcmp(varb,s) == 0),
   tpskip(id,x,y);
   if (feof(id) == 1), err = 1; end;
   [t,s,x,y] = tphead(id);
end

data = tpdata(id,x,y);

fclose(id);
