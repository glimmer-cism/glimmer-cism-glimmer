function [st] = glim2d(want,fname,stem);

% [a] = blah2d('g_land',200e3,'.');
% where a is structure containing 2d arrays of variables
% (eg a.thck is thickness in m)
% 1st param is desired time slice 
% (eg 200 kyr)
% 2nd param is path of results directory relative to current dir
% (eg '.' if blah.2df is in current dir)

% 1 uflx flux in x (m^2/yr)
% 2 vflx flux in y (m^2/yr)
% 3 diffu apparent diffusivity (m^2/yr)
% 4 btrc basal slip coef m/(Pa yr)
% 5 ubas basal slip velocity in x (m/yr)
% 6 vbas basal slip velocity in x (m/yr)
% 7 thck thickness in m
% 8 usrf ice upper surface elevation in m
% 9 lsrf ice lower surface elevation in m
% 10 topg bedrock topography in m
% 11 acab accumulation-ablation rate m/yr
% 12 bmlt basal melt rate in m
% 13 bwat basal water depth m
% 14 artm annual mean air temperature deg. C
% 15 temp(upn) basal ice temperature deg. C (uncorrected)
% 16 arng annual air temperature range deg. C
% 17 prcp precipitation in m/yr
% 18 ablt ablation rate in m/yr
% 19 dusrfdtm rate of upper ice surface elevation change m/yr

% open the file

id = blahopen([fname,'/',stem,'.gl2']);
 
% for staggered grid
% read numbers of dimensions; staggered grid code (0 staggered);
% x/y size for staggered 
 
[dims,stag,ewnv,nsnv,upn] = blahsize(id,2);

uq = zeros(ewnv,nsnv); vq = uq; df = uq; bx = uq; bu = uq; bv = uq;

% read x coords for staggered grid 

[x] = blahcord(id,ewnv);
[y] = blahcord(id,nsnv);

[xs,ys] = meshgrid(y,x); clear x y;

% for normal grid

[dims,norm,ewn,nsn,upn] = blahsize(id,2);

th = zeros(ewn,nsn); up = th; lw = th; to = th; ac = th; bm = th; 
bw = th; at = th; bt = th; ar = th; pc = th; ab = th; ds = th;

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
      case 1; uq = in; disp(['got uq - scale ' num2str(scal)])
      case 2; vq = in; disp(['got vq - scale ' num2str(scal)])
      case 3; df = in; disp(['got df - scale ' num2str(scal)])
      case 4; bx = in; disp(['got bx - scale ' num2str(scal)])
      case 5; bu = in; disp(['got bu - scale ' num2str(scal)])
      case 6; bv = in; disp(['got bv - scale ' num2str(scal)])
    end

  else,

    switch varb
      case 1; th = in; disp(['got th - scale ' num2str(scal)])
      case 2; up = in; disp(['got up - scale ' num2str(scal)])
      case 3; lw = in; disp(['got lw - scale ' num2str(scal)])
      case 4; to = in; disp(['got to - scale ' num2str(scal)])
      case 5; ac = in; disp(['got ac - scale ' num2str(scal)])
      case 6; bm = in; disp(['got bm - scale ' num2str(scal)])
      case 7; bw = in; disp(['got bw - scale ' num2str(scal)])
      case 8; at = in; disp(['got at - scale ' num2str(scal)])
      case 9; bt = in; disp(['got bt - scale ' num2str(scal)])
      case 10; ar = in; disp(['got ar - scale ' num2str(scal)])
      case 11; pc = in; disp(['got pc - scale ' num2str(scal)])
      case 12; ab = in; disp(['got ab - scale ' num2str(scal)])
      case 13; ds = in; disp(['got ds - scale ' num2str(scal)])
    end

  end

  [time,varb,flag,scal] = blahhead(id);

end

fclose(id);

st.xs = xs;
st.ys = ys;
st.uq = uq;
st.vq = vq;
st.df = df;
st.bx = bx;
st.bu = bu;
st.bv = bv;

clear xs ys uq vq df bx bu bv

st.xn = xn;
st.yn = yn;
st.th = th;
st.up = up;
st.lw = lw;
st.to = to;
st.ac = ac;
st.bm = bm;
st.bw = bw;
st.at = at;
st.bt = bt;
st.ar = ar;
st.pc = pc;
st.ab = ab;
st.ds = ds;

clear xn yn th up lw to ac bm bw at bt ar pc ab ds

