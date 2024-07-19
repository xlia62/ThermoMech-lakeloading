% Initial model size (m)
xsize0 = 100e3;
ysize0 = 120e3;

% Number of grid-lines
xnum = 201;
ynum = 61;

% Grid type: uniform, nonuniform, boundary, adaptive (ver construct grid.m)
G.grid_type = 'nonuniform';
% Resolution of fine region (m)
bx = 1000;
by = 1000;
%w(x,y) = width of fine uniform grid (m)
wx = 0;   % Uniform grid in x
wy = 31e3;

%f=center of the fine region i.e f=.5 is dead center f=.25 is a quarter
f = 0.5;% 0.457;