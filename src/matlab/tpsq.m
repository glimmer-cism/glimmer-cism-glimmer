function [o] = tpsq(in,dd,nn)

if dd==1
  o = squeeze(in(nn,:,:));
elseif dd == 2
  o = squeeze(in(:,nn,:));   
elseif dd == 3
  o = squeeze(in(:,:,nn));
end