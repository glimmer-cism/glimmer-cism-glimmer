function [data] = tpread(id,ewn,nsn);
fread(id,1,'int');
data = fread(id,[ewn,nsn],'float');
fread(id,1,'int');
data = flipud(rot90(data));
for pt = 1:ewn*nsn
   if data(pt) < -9999.0 & data(pt) > -10000.0; data(pt) = NaN;
   end
end
