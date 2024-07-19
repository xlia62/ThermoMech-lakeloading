% Initial model size (m)
xsize0 = 1500e3;
ysize0 = 520e3;

% Number of grid-lines
xnum = 251;
ynum = 111;

% Grid type: uniform, nonuniform, boundary, adaptive (ver construct grid.m)
G.grid_type = 'uniform';
% Resolution of fine region (m)
bx = 2000;
by = 2000;
%w(x,y) = width of fine uniform grid (m)
wx = 200e3;
wy = 200e3;
%f=center of the fine region i.e f=.5 is dead center f=.25 is a quarter
f = 0.5;% 0.457;