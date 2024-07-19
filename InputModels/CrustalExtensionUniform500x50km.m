% Initial model size (m)
xsize0 = 500e3;
ysize0 = 70e3;

% Number of grid-lines
xnum = 251;
ynum = 71;

% Grid type: uniform, nonuniform, boundary, adaptive (ver construct grid.m)
G.grid_type = 'nonuniform';
% Resolution of fine region (m)
bx = 1000;
by = 500;
%w(x,y) = width of fine uniform grid (m)
wx = 100e3;
wy = 50e3;

%f=center of the fine region i.e f=.5 is dead center f=.25 is a quarter
f = 0.5;% 0.457;