function out = ustag(in)
[n,m]=size(in);
out=zeros(n+1,m+1);
out(2:end-1,2:end-1) = 0.25 * (in(1:end-1,1:end-1)+in(1:end-1,2:end)+in(2:end,1:end-1)+in(2:end,2:end));
