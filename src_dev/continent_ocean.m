function D = continent_ocean(x,h,b,sea_lvl)

n = length(h);
jnode = 1:n;

% Locate areas above sea level
k = h >= sea_lvl;

% Identify how many segments we need to break topography into
nbrk = sum(abs(diff(k)));

% Generate segment array, first row contains number of points in segment
% second row contains type id: 1 == continent, 0 == ocean
seg = zeros(2,nbrk+1);

j = 1;
seg(1,j) = 1;
seg(2,j) = 1;

for i=2:n
    if seg(2,j) == k(i)
        seg(1,j) = seg(1,j) + 1;
    else
        j = j + 1;
        seg(1,j) = 1;
        seg(2,j) = k(i);
    end
end

D.nseg = nbrk+1;
D.h = mat2cell(h,1,seg(1,:));
D.b = mat2cell(b,1,seg(1,:));
D.x = mat2cell(x,1,seg(1,:));
D.j = mat2cell(jnode,1,seg(1,:));
D.seg = seg;