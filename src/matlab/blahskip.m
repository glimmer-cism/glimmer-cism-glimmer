function blahskip(id,flag,stag,nsnv,nsn,ewnv,ewn,upn)

% this functions skips data in the 3d file

if upn == 0
  if flag == stag  
    skip = 4 * (2 + ewnv * nsnv);
  else
    skip = 4 * (2 + ewn * nsn);
  end
else
  if flag == stag  
    skip = 4 * (2 + ewnv * nsnv * upn);
  else
    skip = 4 * (2 + ewn * nsn * upn);
  end
end
 
fseek(id,skip,0);
