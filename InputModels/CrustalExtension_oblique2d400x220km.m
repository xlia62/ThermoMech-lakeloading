% Initial model size (m)
xsize0 = 400e3; % -> original 400e3
ysize0 = 220e3; % -> original 200e3

% Number of grid-lines
xnum = 485;  % 250 fine of 50 km + 350/1.5 = 50e3/200 fine +234  coarse +1
ynum = 365;  % 250 fine of 50 km + 170/1.5 = 50e3/200 fine +114 coarse +1

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
