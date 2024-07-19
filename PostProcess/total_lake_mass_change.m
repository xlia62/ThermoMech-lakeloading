function [t,Mwater] = total_lake_mass_change(inpath,nstep,npx,npy,xlims,ylims)

% Calculate mass change from lake (kg)

% input:
% Inpath = folder where grid_###.mat and markers_###.mat files exist
% nstep  = time step of interest
% npx,npy = number of image pixels in x and y direction
% xlims,ylims = constrain caclulations to x-limits and y-limits of model, in km

% output: 
% time (Myr) and lake mass change (kg)

% required:
% vis_markers

% main code:
%% load data 

yr2sec = 365.25*24*3600;

% Load files
infile = [inpath,'/markers_',num2str(nstep),'.mat'];
chk = whos('-file',infile);
load(infile);

% 
% if ismember('MXM',{chk.name})
%     load(infile,'MX','MY','MI','MTK', 'MRHO');
% else
%     %MXM = zeros(length(MI));
%     error('No MXM variable found')
% end

infile = [inpath,'/grids_',num2str(nstep),'.mat'];
load(infile);

% Replace other material with NaN
j = MI ~= 2;
MI(j) = NaN;
MX(j) = NaN;
MY(j) = NaN;

if ~isempty(xlims)
    i = MX/1000 >= xlims(1) & MX/1000 <= xlims(2);
    ox = xlims(1);
    mx = MX(i) - ox*1000;
    my = MY(i);
    mi = MI(i);
else
    xlims = [G.gridx(1) G.gridx(end)]/1000;
    ox = 0;
    mx = MX;
    my = MY;
    mi = MI;
end
if ~isempty(ylims)
    i = my/1000 >= ylims(1) & my/1000 <= ylims(2);
    oy = ylims(1);
    my = my(i) - oy*1000;
    mx = mx(i);
    mi = mi(i);
else
    ylims = [G.gridy(1) G.gridy(end)]/1000;
    oy = 0;
end

% 
[~,~,mmass] = vis_markers((xlims-ox)*1000,(ylims-oy)*1000,npx,npy,mx,my,mi);

dx = diff(xlims)/npx; % km
dy = diff(ylims)/npy; % km

% Calculate volume of melt
mmass = mmass(:);
Mwater = sum(mmass(~isnan(mmass))*dx*dy)*MRHO(2)*1e6; % to meter and kg

% Corresponding time
t = timesum*1e-6/(yr2sec);
