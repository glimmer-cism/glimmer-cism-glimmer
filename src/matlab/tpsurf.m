function tpsurf(x,z,a,dim,num)

if dim == 3
  surf(x(:,num),tpsq(z,dim,num),tpsq(a,dim,num));
elseif dim ==2
  surf(x(num,:),tpsq(z,dim,num),tpsq(a,dim,num));
end

view(2); axis tight; colorbar