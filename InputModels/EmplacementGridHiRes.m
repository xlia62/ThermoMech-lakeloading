% Initial model size (m)
xsize0 = 100e3;
ysize0 = 140e3;

% Number of grid-lines
xnum = 251;
ynum = 141;

% Grid type: uniform, nonuniform, boundary, adaptive (ver construct grid.m)
G.grid_type = 'nonuniform';
% Resolution of fine region (m)
bx = 500;
by = 500;
%w(x,y) = width of fine uniform grid (m)
wx = 0;   % Uniform grid in x
wy = 45e3;

%f=center of the fine region i.e f=.5 is dead center f=.25 is a quarter
f = 0.5;% 0.457;