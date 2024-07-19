function quad = connect_quad(nx,ny)

nex = nx - 1;
ney = ny - 1;

quad = zeros(nex*ney,4);

quad(1,:) = [1 2 2+nx 1+nx];

for i = 2:nex
    quad(i,:) = quad(1,:) + i - 1;
end

j = nex + 1;
k = j + nex -1;

for i = 1:ney-1
    quad(j:k,:) = quad(1:nex,:) + i*nx;
    j = j + nex;
    k = k + nex;
end