function [lon,lat,gds,bas] = forcebal(a,i)

lat = (mean(a.txy(:,:,i+1),1).*(a.zs(1,:,i+1)-a.zs(end,:,i+1)) - ...
       mean(a.txy(:,:,i),1).*(a.zs(1,:,i)-a.zs(end,:,i))) ./ ...
       (a.ys(end,i+1)-a.ys(end,i));

lon = lat; lon(1) = 0.0; lon(end) = 0.0;

lon(2:end-1) = 2*(mean(a.txx(:,3:end,i+1),1).*(a.zs(1,3:end,i+1)-a.zs(end,3:end,i+1)) - ...
       mean(a.txx(:,1:end-2,i),1).*(a.zs(1,1:end-2,i)-a.zs(end,1:end-2,i))) ./ ...
       (a.xs(end,i)-a.xs(end-2,i));

gds = mean(a.gdx(end,:,i:i+1),3);

bas = mean(a.txz(end,:,i:i+1),3);

xx = [1:size(bas,2)];

plot(xx,-lon,'g',xx,-lat,'y',xx,bas,'m',xx,gds,'r-o',xx,-lon-lat+bas,'k+-')
legend('lon','lat','bas','gds','tot')