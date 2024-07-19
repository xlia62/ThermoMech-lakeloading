% Initial model size (m)
xsize0 = 1000e3;
ysize0 = 700e3;

% Number of grid-lines
xnum = 201;
ynum = 141;

% Grid type: uniform, nonuniform, boundary, adaptive (ver construct grid.m)
G.grid_type = 'uniform';

% Not used, but code needs it (to be fixed)
bx = 2000;
by = 2000;
%w(x,y) = width of fine uniform grid (m)
wx = 300e3;
wy = 100e3;
%f=center of the fine region i.e f=.5 is dead center f=.25 is a quarter
f = 0.5;% 0.457;