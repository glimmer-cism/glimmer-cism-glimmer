function [x] = blahcord(id,sz)

fread(id,1,'int');
[x] = fread(id,[sz],'float');
fread(id,1,'int');

