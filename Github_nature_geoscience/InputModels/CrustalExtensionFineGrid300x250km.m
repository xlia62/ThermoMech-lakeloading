% Initial model size (m)
xsize0 = 300e3; % -> original 400e3
ysize0 = 250e3; % -> original 200e3

% Number of grid-lines
xnum = 418;  % 250 fine + 250/1.5 = 167  coarse +1
ynum = 384;  % 250 fine + 200/1.5 = 133 coarse +1

% Grid type: uniform, nonuniform, boundary, adaptive (ver construct grid.m)
G.grid_type = 'nonuniform';
% Resolution of fine region (m)
bx = 200;
by = 200;
%w(x,y) = width of fine uniform grid (m)
wx = 50e3; 
wy = 50e3; % -> 50e3

%f=center of the fine region i.e f=.5 is dead center f=.25 is a quarter
f = 0.5;% 0.457;
