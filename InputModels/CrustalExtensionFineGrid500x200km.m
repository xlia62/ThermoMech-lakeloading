% Initial model size (m)
xsize0 = 400e3; % -> original 400e3
ysize0 = 200e3; % -> original 200e3

% Number of grid-lines
xnum = 421;
ynum = 271;

% Grid type: uniform, nonuniform, boundary, adaptive (ver construct grid.m)
G.grid_type = 'nonuniform';
% Resolution of fine region (m)
bx = 250;
by = 250;
%w(x,y) = width of fine uniform grid (m)
wx = 50e3; 
wy = 50e3; % -> 50e3

%f=center of the fine region i.e f=.5 is dead center f=.25 is a quarter
f = 0.5;% 0.457;