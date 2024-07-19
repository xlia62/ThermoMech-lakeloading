% Initial model size (m)
%xsize0 = 2700e3 + 150e3;  % in SUBDUCT20, widened shelf by 150 km
%xsize0 = 3000e3; % try for subduct 24 (might have been probelem in subduct 23)
%ysize0 = 1000e3;
%ysize0 = 1000e3; % also for subduct 24
xsize0 = 2700e3 + 150e3;  % in SUBDUCT20, widened shelf by 150 km
ysize0 = 1500e3;

% Number of grid-lines
xnum = 251;
ynum = 426;

% Grid type: uniform, nonuniform, boundary, adaptive (ver construct grid.m)
G.grid_type = 'nonuniform';
% Resolution of fine region (m)
bx = 2000;
by = 2000;
%w(x,y) = width of fine uniform grid (m)
wx = 300e3;
wy = 100e3;
%f=center of the fine region i.e f=.5 is dead center f=.25 is a quarter
f = 0.52;% 0.457;