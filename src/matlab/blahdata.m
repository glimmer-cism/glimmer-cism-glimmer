function [var] = blahdata(id,flag,stag,nsnv,nsn,ewnv,ewn,upn)

% this functions reads data from 3d file

% read the variable according whether it is staggered grid or not

fread(id,1,'int');

if upn == 0
  if flag == stag  
    var= fread(id,[nsnv,ewnv],'float');
  else
    var = fread(id,[nsn,ewn],'float');
  end
  var = var';
else
  if flag == stag  
    var = fread(id,[upn*ewnv*nsnv],'float'); 
    var = reshape(var,upn,nsnv,ewnv);
  else
    var = fread(id,[upn*ewn*nsn],'float'); 
    var = reshape(var,upn,nsn,ewn);
  end
end

fread(id,1,'int');


