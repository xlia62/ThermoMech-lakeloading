% Initial model size (m)
xsize0 = 3000e3;
ysize0 = 750e3;

% Number of grid-lines
xnum = 331;
ynum = 171;

% Grid type: uniform, nonuniform, boundary, adaptive (ver construct grid.m)
G.grid_type = 'nonuniform';
% Resolution of fine region (m)
bx = 2000;
by = 2000;
%w(x,y) = width of fine uniform grid (m)
wx = 250e3;
wy = 250e3;
%f=center of the fine region i.e f=.5 is dead center f=.25 is a quarter
f = 0.65; %0.5 