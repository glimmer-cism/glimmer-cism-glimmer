function id = blahopen(filename);

% this function opens a blah file and checks to see
% if it is in a readable format, if not it opens
% it in a different binary format

id = fopen(filename,'r','ieee-be');
a = fread(id,1,'int');

if (a > 100)
  fclose(id);
  id = fopen(filename,'r','ieee-le');
else
  frewind(id);
end

