function T=plot_t_meltV_reference(inpath,nbegin, nfreq, nend,xlims,ylims, melt_i)
% this is function plot the melting history_for_refernece(km^3/km)

% input:
% Inpath = folder where grid_###.mat and markers_###.mat files exist
% nfreq  = frequency of time stpes when files are saved
% nbegin = start number of time steps 
% nend = end number of time steps
% npx,npy = number of image pixels in x and y direction
% melt choose to calculate melt or extract 1:MXM or 2:MEXTC

% xlims,ylims = constrain caclulations to x-limits and y-limits of model, in km
% rho_solid = density of the solid phase kg/m^3
% rho_melt_i = density of the melt phase kg/m^3
% 9 = Continental Lithospheric mantle -- dry mantle 
% MRHO(9,1) = 3300;            % standard density, kg/m^3
% MRHO(9,4) = 2900;            % melt density, kg/m^3 
% 7 = Lower continental crust (diorite)
% MRHO(7,1) = 2900;             % standard density, kg/m^3
% MRHO(7,4) = 2400;             % melt density, kg/m^3
% Here is the default value of density
rho_solid = 3300;
rho_melt = 2900;
marker_per_grid = 36;

% output: 
% image of melting history(km^3/km)

% required:
% totalmelt_lake function

% main code:
%% load data and define indexes
% if nfreq > 1
%     k = [1 nfreq:nfreq:endstep];
% else
%     k = [nfreq:nfreq:endstep];
% end
k = [nbegin:nfreq:nend];

% Load files
infile = [inpath,'/markers_',num2str(nend),'.mat'];
load(infile);
infile = [inpath,'/grids_',num2str(nend),'.mat'];
load(infile);
figname = split(inpath, '/'); %
%refername = replace(char(figname(end)), 'drop', 'nochange'); 


% Replace Air or Water with NaN
i = MI == 1 | MI == 2;
MGII(i) = NaN;

% select based on x and y limits
if ~isempty(xlims)
    i = MX/1000 >= xlims(1) & MX/1000 <= xlims(2);
    ox = xlims(1);
    mx = MX(i) - ox*1000;
    my = MY(i);
% find out resolution of marker
    grid_x_num = G.xnum * (xlims(2)-xlims(1))*1e3/(max(MX) - min(MX)); % number of x grids
    npx = round(grid_x_num) *marker_per_grid; % 4 markers in each grid   
else
    xlims = [G.gridx(1) G.gridx(end)]/1000;
    ox = 0;
    mx = MX;
    my = MY;
    npx = G.xnum*marker_per_grid;
end

if ~isempty(ylims)
    i = my/1000 >= ylims(1) & my/1000 <= ylims(2);
    oy = ylims(1);
    my = my(i) - oy*1000;
    mx = mx(i);
% find out resolution of marker
    grid_y_num = G.ynum * (ylims(2)-ylims(1))*1e3/(max(MY) - min(MY)); % number of x grids
    npy = round(grid_y_num) *marker_per_grid; % 4 markers in each grid   
else
    ylims = [G.gridy(1) G.gridy(end)]/1000;
    oy = 0;
    my = my;
    npy = G.ynum*marker_per_grid;
end


time = zeros(length(k), 1);
melt = zeros(length(k), 1);
melt_drop = zeros(length(k), 1);
water = zeros(length(k), 1);
melt_fraction = zeros(length(k), 1);
% loop to 
count = 1;
for i = k(1:end)
    if melt_i == 1
        [t_i,Vmelt_i_drop] = totalmelt_lake(inpath,i,npx,npy,xlims,ylims);
        % please note the we add a reference of melting 
        %[t_i,Vmelt_i] = totalmelt_lake(['output/AfricaModels2022/Lake_model_extraction/', refername],i,npx,npy,xlims,ylims);
        type = 'melting';
    elseif melt_i == 2
        [t_i,Vmelt_i_drop] = totalmeltextracted(inpath,i,npx,npy,xlims,ylims, 3300, 2900);
        % please note the we add a reference of melting 
        %[t_i,Vmelt_i] = totalmeltextracted(['output/AfricaModels2022/Lake_model_extraction/', refername],i,npx,npy,xlims,ylims, 3300, 2900);
        type = 'melting_extract';
    else
        display('Wrong input of melting type.')
    end
    [t_i,Mwater_i] = total_lake_mass_change(inpath,i,npx,npy,xlims,ylims);
    time(count) = t_i;
    % find the max 10 particles and average mf 
    mf_infile = [inpath,'/markers_',num2str(i),'.mat'];
    mf = load(mf_infile);
    A = sort(mf.MXM);
    mf_i = sum(A([end-5: end]))/6; 
    melt_drop(count) = Vmelt_i_drop;
    %melt(count) = Vmelt_i;
    water(count) = Mwater_i;
    melt_fraction(count) = mf_i;
    count = count +1; 
end

time_ky = time*1e3;
melt_v = melt_drop; %km^3/km
water_mass_kg = water; %kg
max_mf = melt_fraction;  % max accumulated melt fraction
T = table(time_ky, melt_v,  water_mass_kg, max_mf);
output_table = ['output/AfricaModels2022/Table/',char(figname(end)),'_',...
    num2str(type),'_history_',num2str(nbegin), '-', num2str(nend), ...
    '_dt_', num2str(nfreq), '_nomark_', num2str(marker_per_grid), ...
    '_x_', num2str(xlims(1)), '-',  num2str(xlims(2)),...
    '_y_', num2str(ylims(1)), '-',  num2str(ylims(2))];
writetable(T, [output_table, '.csv']); 


end
