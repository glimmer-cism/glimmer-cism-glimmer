function [data] = tpread(id,ewn,nsn);
skip = 4 * (2 + ewn * nsn);
fseek(id,skip,0);
