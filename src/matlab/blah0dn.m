function [st] = blah0dn(fname);

% [st] = blah0dn('.');
% [st] = blah0dn('.');

nv = 14; np = 4;

id = blahopen([fname,'/blah.0df']);

fread(id,1,'int'); dims = fread(id,1,'int'); fread(id,1,'int');

fread(id,1,'int'); ewv = fread(id,dims,'float'); fread(id,1,'int');
fread(id,1,'int'); nsv = fread(id,dims,'float'); fread(id,1,'int');
fread(id,1,'int'); ew = fread(id,dims,'float'); fread(id,1,'int');
fread(id,1,'int'); ns = fread(id,dims,'float'); fread(id,1,'int');

fread(id,1,'int');
[in c] = fread(id,np+nv*dims,'float');  
fread(id,1,'int');  

i = 1;

while c == np+nv*dims
 
  tm(i) = in(1);
  vo(i) = in(2);
  ar(i) = in(3);
  mt(i) = in(4); 

  s = np+1; e = np+dims;

  uq(1:dims,i) = in(s:e); s = s + dims; e = e + dims;
  vq(1:dims,i) = in(s:e); s = s + dims; e = e + dims;
  bu(1:dims,i) = in(s:e); s = s + dims; e = e + dims;
  bv(1:dims,i) = in(s:e); s = s + dims; e = e + dims;
  bx(1:dims,i) = in(s:e); s = s + dims; e = e + dims;
  th(1:dims,i) = in(s:e); s = s + dims; e = e + dims;
  up(1:dims,i) = in(s:e); s = s + dims; e = e + dims;
  lw(1:dims,i) = in(s:e); s = s + dims; e = e + dims;
  to(1:dims,i) = in(s:e); s = s + dims; e = e + dims;
  ac(1:dims,i) = in(s:e); s = s + dims; e = e + dims;
  bm(1:dims,i) = in(s:e); s = s + dims; e = e + dims;    
  bw(1:dims,i) = in(s:e); s = s + dims; e = e + dims;
  at(1:dims,i) = in(s:e); s = s + dims; e = e + dims;
  bt(1:dims,i) = in(s:e);

  i = i + 1;

  fread(id,1,'int');
  [in c] = fread(id,np+nv*dims,'float');  
  fread(id,1,'int');  

end

fclose(id);

st.tm = tm;
st.vo = vo;
st.ar = ar;
st.mt = mt;
st.uq = uq;
st.vq = vq;
st.bu = bu;
st.bv = bv;
st.bx = bx;
st.th = th;
st.up = up;
st.lw = lw;
st.to = to;
st.ac = ac;
st.bm = bm;
st.bw = bw;
st.at = at;
st.bt = bt;

st.ewv = ewv;
st.ew = ew;
st.nsv = nsv;
st.ns = ns;

clear tm ewv nsv ew ns vo ar mt uq vq bu bv bx th up lw to ac bm bw at bt






