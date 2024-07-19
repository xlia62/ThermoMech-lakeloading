% Initial model size (m)
xsize0 = 4000e3;
ysize0 = 800e3;

% Number of grid-lines
xnum = 401;
ynum = 151;

% Grid type: uniform, nonuniform, boundary, adaptive (ver construct grid.m)
G.grid_type = 'gerya_subduction';
% xpos = x position of the start of the fine region (for gerya_subduction)
G.xpos = 2000e3;

% Resolution of fine region (m)
bx = 2000;
by = 2000;
%w(x,y) = width of fine uniform grid (m)
wx = 250e3;
wy = 120e3;
%f=center of the fine region i.e f=.5 is dead center f=.25 is a quarter
f = 0.5;% 0.457;


% Note: move this high resolution area
% To view the grids
% load('grids_##.mat')
% [pgx,pgy] = plotgrid(G.gridx/1000, G.gridy/1000)
% plot(pgx(:,1),pgx(:,2))
% plot(pgy(:,1),pgy(:,2))