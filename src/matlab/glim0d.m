function [tm,vo,ar,mt,uq,vq,th,up,lw,to,ac,bm,bw,at,bt] = glim0d(fname,stem);

id = fopen([fname,'/',stem,'.gl0']);
 
fread(id,1,'int'); dims = fread(id,1,'int'); fread(id,1,'int');

fread(id,1,'int'); ewv = fread(id,dims,'float'); fread(id,1,'int');
fread(id,1,'int'); nsv = fread(id,dims,'float'); fread(id,1,'int');
fread(id,1,'int'); ew = fread(id,dims,'float'); fread(id,1,'int');
fread(id,1,'int'); ns = fread(id,dims,'float'); fread(id,1,'int');

fread(id,1,'int');
[in c] = fread(id,4+11*dims,'float');  
fread(id,1,'int');  

i = 1;

while c == 4+11*dims
 
  tm(i) = in(1);
  vo(i) = in(2);
  ar(i) = in(3);
  mt(i) = in(4); 

  uq(1:dims,i) = in(5:4+dims); 
  vq(1:dims,i) = in(5+dims:4+2*dims);
  th(1:dims,i) = in(5+2*dims:4+3*dims);
  up(1:dims,i) = in(5+3*dims:4+4*dims); 
  lw(1:dims,i) = in(5+4*dims:4+5*dims); 
  to(1:dims,i) = in(5+5*dims:4+6*dims); 
  ac(1:dims,i) = in(5+6*dims:4+7*dims);
  bm(1:dims,i) = in(5+7*dims:4+8*dims); 
  bw(1:dims,i) = in(5+8*dims:4+9*dims); 
  at(1:dims,i) = in(5+9*dims:4+10*dims); 
  bt(1:dims,i) = in(5+10*dims:4+11*dims);

  i = i + 1;

  fread(id,1,'int');
  [in c] = fread(id,4+11*dims,'float');  
  fread(id,1,'int');  

end

fclose(id);




