function [time,str,ewn,nsn] = tpread(id);
fread(id,1,'int');
time = fread(id,1,'float');
str = setstr(rot90(fread(id,4,'char')));
ewn = fread(id,1,'int');
nsn = fread(id,1,'int');
dxy = fread(id,1,'float');
fread(id,1,'int');
