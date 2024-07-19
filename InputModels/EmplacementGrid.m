% Initial model size (m)
xsize0 = 100e3;
ysize0 = 200e3;

% Number of grid-lines
xnum = 321; % was 201
ynum = 101; % was 63

% Grid type: uniform, nonuniform, boundary, adaptive (ver construct grid.m)
G.grid_type = 'nonuniform';
% Resolution of fine region (m)
bx = 1000;
by = 1000;
%w(x,y) = width of fine uniform grid (m)
wx = 0;   % Uniform grid in x
wy = 44e3; % was 31e3 added 13e3 to account for increased sticky layer

%f=center of the fine region i.e f=.5 is dead center f=.25 is a quarter
f = 0.5;% 0.457;