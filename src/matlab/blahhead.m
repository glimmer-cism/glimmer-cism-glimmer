function [time,varb,flag,scal] = blahhead(id)

% this functions reads data header from 3d file

% read header which is time stamp; variable code; grid-type code; scaling
 
fread(id,1,'int');
time = fread(id,1,'float'); varb = fread(id,1,'int'); flag = fread(id,1,'int'); scal = fread(id,1,'float'); 
fread(id,1,'int');
