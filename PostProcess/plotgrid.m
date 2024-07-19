function [pgx,pgy]=plotgrid(gx,gy)

nx = length(gx);
ny = length(gy);

gxx = [gx(:) gx(:) NaN*gx(:)]';
gxy = [gy(1) gy(end) NaN]';
gyy = [gy(:) gy(:) NaN*gy(:)]';
gyx = [gx(1) gx(end) NaN]';

gxx= gxx(:);
gxy = repmat(gxy(:),nx,1);
gyy = gyy(:);
gyx = repmat(gyx(:),ny,1);

pgx = [gxx gxy];
pgy = [gyx gyy];

if nargout == 0
    figure
    plot(pgx(:,1),pgx(:,2),'k')
    hold on
    plot(pgy(:,1),pgy(:,2),'k')
    axis image
end