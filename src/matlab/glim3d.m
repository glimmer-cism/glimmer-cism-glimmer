function [st] = glim3d(want,fname,stem);

% [a] = blah3d(200e3,'.');
% [a] = blah3d(200e3,'.');

% open the file

id = blahopen([fname,'/',stem,'.gl3']);
 
% for staggered grid
% read numbers of dimensions; staggered grid code (0 staggered);
% x/y size for staggered 
 
[dims,stag,ewnv,nsnv,upn] = blahsize(id,3);

zs = zeros(upn,ewnv,nsnv); u = zs; v = zs;
 
% read x coords for staggered grid 

[x] = blahcord(id,ewnv);
[y] = blahcord(id,nsnv);

[xs,ys] = meshgrid(y,x); clear x y;

% for normal grid

[dims,norm,ewn,nsn,upn] = blahsize(id,3);

zn = zeros(upn,ewn,nsn); w = zn; wg = zn; a = zn; t = zn; 

% read x/y coords for normal grid 

[x] = blahcord(id,ewn);
[y] = blahcord(id,nsn);

[xn,yn] = meshgrid(y,x); clear x y;
 
% now on to variables - staggered grid first

% header and data for u flux

[time,varb,flag,scal] = blahhead(id);

while want ~= time,

  blahskip(id,flag,stag,nsnv,nsn,ewnv,ewn,upn);
  [time,varb,flag,scal] = blahhead(id);   

end

while want == time,

  in = blahdata(id,flag,stag,nsnv,nsn,ewnv,ewn,upn);

  if flag == stag,

    switch varb
      case 1; zs = in; disp(['got zs - scale ' num2str(scal)])
      case 2; u = in; disp(['got u - scale ' num2str(scal)])
      case 3; v = in; disp(['got v - scale ' num2str(scal)])
      case 4; evs = in; disp(['got evs - scale ' num2str(scal)])
      case 5; tau = in; disp(['got atu - scale ' num2str(scal)])
      case 6; txz = in; disp(['got txz - scale ' num2str(scal)])
      case 7; tyz = in; disp(['got tyz - scale ' num2str(scal)])
      case 8; txy = in; disp(['got txy - scale ' num2str(scal)])
      case 9; txx = in; disp(['got txx - scale ' num2str(scal)])
      case 10; tyy = in; disp(['got tyy - scale ' num2str(scal)])
      case 11; gdx = in; disp(['got gdx - scale ' num2str(scal)])
      case 12; gdy = in; disp(['got gdy - scale ' num2str(scal)])

    end

  else,

    switch varb
      case 1; zn = in; disp(['got zn - scale ' num2str(scal)])
      case 2; w = in; disp(['got w - scale ' num2str(scal)])
      case 3; wg = in; disp(['got wg - scale ' num2str(scal)])
      case 4; a = in; disp(['got a - scale ' num2str(scal)])
      case 5; t = in; disp(['got t - scale ' num2str(scal)])
    end

  end

  [time,varb,flag,scal] = blahhead(id);

end

fclose(id);

% tt = t;
% for up = 1:upn; tt(up,:,:) = t(up,:,:) + 9.81 * 910 * 9.8e-8 *(zn(1,:,:) - zn(up,:,:)); end

st.xs = xs';
st.ys = ys';
st.zs = zs;
st.xn = xn';
st.yn = yn';
st.zn = zn;
st.u = u;
st.v = v;
st.evs = evs;
st.tau = tau;
st.txz = txz;
st.tyz = tyz;
st.txy = txy;
st.txx = txx;
st.tyy = tyy;
st.gdx = gdx;
st.gdy = gdy;
st.w = w;
st.wg = wg;
st.a = a;
st.t = t;

clear xs ys zs xn yn zn a t u v w wg


