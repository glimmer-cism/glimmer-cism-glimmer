function [dims,flag,nsnx,ewnx,upnx] = blahsize(id,which)

% this function reads the size info at the start of a 3d file

% read numbers of dimensions (3); grid code (0 staggered, 1 normal);
% x size for grid; y size for grid

fread(id,1,'int');
dims = fread(id,1,'int'); flag = fread(id,1,'int');
nsnx = fread(id,1,'int'); ewnx = fread(id,1,'int');
if which==3; upnx = fread(id,1,'int'); else; upnx = 0; end 
fread(id,1,'int');
