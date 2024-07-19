% Initial model size (m)
xsize0 = 300e3; % -> original 400e3
ysize0 = 250e3; % -> original 200e3

% Number of grid-lines
xnum = 501;  % 50e3/150 fine + 250/1.5 = 333 fine +167  coarse +1
ynum = 467;  % 50e3/150 fine + 200/1.5 = 333 fine +133 coarse +1

% Grid type: uniform, nonuniform, boundary, adaptive (ver construct grid.m)
G.grid_type = 'nonuniform';
% Resolution of fine region (m)
bx = 150;
by = 150;
%w(x,y) = width of fine uniform grid (m)
wx = 50e3; 
wy = 50e3; % -> 50e3

%f=center of the fine region i.e f=.5 is dead center f=.25 is a quarter
f = 0.5;% 0.457;
