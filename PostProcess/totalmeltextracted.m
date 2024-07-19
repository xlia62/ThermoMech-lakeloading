function [t,Vmelt] = totalmeltextracted(inpath,nstep,npx,npy,xlims,ylims,rho_solid,rho_melt)

% Calculate cumulative volume of melt extracted (km^3/km)
%
% [t,Vmelt] = totalmeltextracted(inpath,nstep,npx,npy,xlims,ylims,rho_solid,rho_melt)
%
% Inpath = folder where grid_###.mat and markers_###.mat files exist
% nstep  = time step of interest
% npx,npy = number of image pixels in x and y direction
% xlims,ylims = constrain caclulations to x-limits and y-limits of model, in km
% rho_solid = density of the solid phase kg/m^3
% rho_melt = density of the melt phase kg/m^3
%
% V1.0 Robert Moucha, Syracuse University
%

yr2sec = 365.25*24*3600;

% Load files
infile = [inpath,'/markers_',num2str(nstep),'.mat'];
chk = whos('-file',infile);

if ismember('MEXTC',{chk.name})
    load(infile,'MX','MY','MI','MTK','MEXTC', 'MEXT');
else
    error('No MEXTC variable found')
end

infile = [inpath,'/grids_',num2str(nstep),'.mat'];
load(infile);

if ~isempty(xlims)
    i = MX/1000 >= xlims(1) & MX/1000 <= xlims(2);
    ox = xlims(1);
    mx = MX(i) - ox*1000;
    my = MY(i);
    mtemp = MTK(i);
    mextc = MEXTC(i);
    mext = MEXT(i);
    mi = MI(i);
else
    xlims = [G.gridx(1) G.gridx(end)]/1000;
    ox = 0;
    mx = MX;
    my = MY;
    mtemp = MTK;
    mextc = MEXTC;
    mext = MEXT;
    mi = MI;
end
if ~isempty(ylims)
    i = my/1000 >= ylims(1) & my/1000 <= ylims(2);
    oy = ylims(1);
    my = my(i) - oy*1000;
    mx = mx(i);
    mtemp = mtemp(i);
    mextc = mextc(i);
    mext = mext(i);
    mi = mi(i);
else
    ylims = [G.gridy(1) G.gridy(end)]/1000;
    oy = 0;
end

% Compute cumulative fraction of melt extrected
[~,~,mmass] = vis_markers((xlims-ox)*1000,(ylims-oy)*1000,npx,npy,mx,my,mextc);

dx = diff(xlims)/npx; % in km 
dy = diff(ylims)/npy; % in km 

% Calculate volume of melt
mmass = mmass(:);
Vmelt = sum(mmass(~isnan(mmass))*dx*dy)*rho_melt/rho_solid; % in km^2 

% Corresponding time
t = timesum*1e-6/(yr2sec);
