% % Initial model size (m)
% xsize0 = 6000e3;
% ysize0 = 1500e3;
% 
% % Number of grid-lines
% %xnum = 291;
% xnum = 521;
% ynum = 131;
% % xnum = 781;
% % ynum = 256;

% Initial model size (m)
xsize0 = 3000e3;
ysize0 = 1000e3;

% Number of grid-lines
%xnum = 291;
xnum = 301;
ynum = 101;
% xnum = 781;
% ynum = 256;

% Grid type: uniform, nonuniform, boundary, adaptive (ver construct grid.m)
G.grid_type = 'gerya_subduction';
% Resolution of fine region (m)
bx = 2000;
by = 2000;
%w(x,y) = width of fine uniform grid (m)
wx = 200e3;
wy = 62e3;
%f=center of the fine region i.e f=.5 is dead center f=.25 is a quarter
f = 0.5;% 0.457;
% % xpos = x position of the start of the fine region (for gerya_subduction)
% xpos = 200e3;
% Initial trench position 
init_trenchpos = xsize0*2/3; % km
G.xpos = init_trenchpos;
% init_trenchpos = xsize0*2/3; % km

% Note: move this high resolution area
% To view the grids
% load('grids_##.mat')
% [pgx,pgy] = plotgrid(G.gridx/1000, G.gridy/1000)
% plot(pgx(:,1),pgx(:,2))
% plot(pgy(:,1),pgy(:,2))