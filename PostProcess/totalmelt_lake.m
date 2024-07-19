function [t,Vmelt] = totalmeltextracted_lake(inpath,nstep,npx,npy,xlims,ylims)

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
% 
% 9 = Continental Lithospheric mantle -- dry mantle 
% MRHO(9,1) = 3300;            % standard density, kg/m^3
% MRHO(9,4) = 2900;            % melt density, kg/m^3 
% 7 = Lower continental crust (diorite)
% MRHO(7,1) = 2900;             % standard density, kg/m^3
% MRHO(7,4) = 2400;             % melt density, kg/m^3
rho_solid = 3300;
rho_melt = 2900;


yr2sec = 365.25*24*3600;

% Load files
infile = [inpath,'/markers_',num2str(nstep),'.mat'];
chk = whos('-file',infile);

if ismember('MXM',{chk.name})
    load(infile,'MX','MY','MI','MTK','MXM');
else
    %MXM = zeros(length(MI));
    error('No MXM variable found')
end

infile = [inpath,'/grids_',num2str(nstep),'.mat'];
load(infile);

if ~isempty(xlims)
    i = MX/1000 >= xlims(1) & MX/1000 <= xlims(2);
    ox = xlims(1);
    mx = MX(i) - ox*1000;
    my = MY(i);
    mtemp = MTK(i);
    mxm = MXM(i);
    mi = MI(i);
else
    xlims = [G.gridx(1) G.gridx(end)]/1000;
    ox = 0;
    mx = MX;
    my = MY;
    mtemp = MTK;
    mxm = MXM;
    mi = MI;
end
if ~isempty(ylims)
    i = my/1000 >= ylims(1) & my/1000 <= ylims(2);
    oy = ylims(1);
    my = my(i) - oy*1000;
    mx = mx(i);
    mtemp = mtemp(i);
    mxm = mxm(i);
    mi = mi(i);
else
    ylims = [G.gridy(1) G.gridy(end)]/1000;
    oy = 0;
end

% Compute cumulative fraction of melt extrected
[~,~,mmass] = vis_markers((xlims-ox)*1000,(ylims-oy)*1000,npx,npy,mx,my,mxm);

dx = diff(xlims)/npx; % km
dy = diff(ylims)/npy; % km

% Calculate volume of melt
mmass = mmass(:);
Vmelt = sum(mmass(~isnan(mmass))*dx*dy); % km^2

% Corresponding time
t = timesum*1e-6/(yr2sec);
