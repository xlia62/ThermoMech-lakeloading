function [x,y] = quassidistribute(width,height,nx,ny,weight)

totalArea = width*height;
pointArea = totalArea/(nx*ny);

L = sqrt(pointArea);

vi = linspace(0,width,2*nx+1);
vj = linspace(0,height,2*ny+1);
vi = vi(2:2:end-1);
vj = vj(2:2:end-1);

dx = vi(1);

% Reduce dy if only one value (change weight in y)
dy = vj(1);
if ny == 1
    dy = 0.9*dy;
end

x = zeros(nx*ny,1); y = x;

c = 0;
for i = 1:nx
    for j = 1:ny
        c = c+1;
        x(c) = vi(i) + (dx*(rand(1))-dx/2)*weight;
        y(c) = vj(j) + (dy*(rand(1))-dy/2)*weight;
    end
end