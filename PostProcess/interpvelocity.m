function [Vx,Vy]=interpvelocity(X,Y,G)

% Interpolating vx, vy for the basic grid

vxb=zeros(G.ynum,G.xnum);
vyb=zeros(G.ynum,G.xnum);

for j=1:G.xnum
    for i=1:G.ynum
        vxb(i,j)=(G.vx1(i,j)+G.vx1(i+1,j))/2;
        vyb(i,j)=(G.vy1(i,j)+G.vy1(i,j+1))/2;
    end
end

Vx = interp2(G.gridx,G.gridy,vxb,X,Y);
Vy = interp2(G.gridx,G.gridy,vyb,X,Y);


