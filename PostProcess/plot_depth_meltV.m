function T=plot_depth_meltV(inpath,step, nfreq, xlims,ylims)
% this is function plotting the extracted melting volume vs depth (km^3/km)

% input:
% Inpath = folder where grid_###.mat and markers_###.mat files exist
% step: the step in which melting volume should be calculated 
% nfreq: the numbers of depth interval of calculation 
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
factor = 1.01;
% output: 
% image of melting history(km^3/km)

% required:
% totalmelt_lake function

% main code:
%% load data and define indexes

depth = linspace(20, ylims(2), nfreq+1); % 20 km is the thickness of the sticky air 
interval =  (ylims(2)-20)/nfreq; %depth interval
% Load files
infile = [inpath,'/markers_',num2str(step),'.mat'];
load(infile);
infile = [inpath,'/grids_',num2str(step),'.mat'];
load(infile);
figname = split(inpath, '/'); 
refername = replace(char(figname(end)), 'drop', 'nochange'); 


melt = zeros(length(depth), 1);
melt_drop = zeros(length(depth), 1); %(km^3/km)
depth_z = zeros(length(depth), 1); % depth in km 
melt_percentage = zeros(length(depth), 1); %(km^3/km)
melt_diff = zeros(length(depth), 1); %(km^3/km)
count = 1;

for d = depth 
    % Replace Air or Water with NaN
    i = MI == 1 | MI == 2;
    MGII(i) = NaN;
    
    depthlims = [d, d+interval];
    
    % select based on x and y limits
    if ~isempty(xlims)
        i = MX/1000 >= xlims(1) & MX/1000 <= xlims(2);
        ox = xlims(1);
        mx = MX(i) - ox*1000;
        my = MY(i);
    % find out resolution of marker, npx, npy
        grid_x_num = G.xnum * (xlims(2)-xlims(1))*1e3/(max(MX) - min(MX)); % number of x grids
        npx = round(grid_x_num) *marker_per_grid; % using markers in each grid   
    else
        xlims = [G.gridx(1) G.gridx(end)]/1000;
        ox = 0;
        mx = MX;
        my = MY;
        npx = G.xnum*marker_per_grid;
    end
    % now calculate each depth interval
    if ~isempty(depthlims)
        i = my/1000 >= depthlims(1) & my/1000 <= depthlims(2);
        oy = depthlims(1);
        my = my(i) - oy*1000;
        mx = mx(i);
    % find out resolution of marker
        grid_y_num = G.ynum * (depthlims(2)-depthlims(1))*1e3/(max(MY) - min(MY)); % number of x grids
        npy = round(grid_y_num) *marker_per_grid; % 4 markers in each grid   
    else
        ylims = [G.gridy(1) G.gridy(end)]/1000;
        oy = 0;
        my = my;
        npy = G.ynum*marker_per_grid;
    end
    
    
    [t_i,Vmelt_i_drop] = totalmeltextracted(inpath,step,npx,npy,xlims,depthlims, rho_solid, rho_melt); % km^3/km
    % please note the we add a reference of melting 
    [t_i,Vmelt_i] = totalmeltextracted(['output/AfricaModels2022/Lake_model_extraction/', refername],step,npx,npy,xlims,depthlims, rho_solid, rho_melt);
    pressure_diff = Vmelt_i_drop- Vmelt_i;
    if pressure_diff <= 0
        pressure_pre = 0;
    else
        pressure_pre = 100* pressure_diff/Vmelt_i;
    end 

% write all data
    melt_drop(count) = Vmelt_i_drop*factor;
    melt(count) = Vmelt_i;
    depth_z(count) = d- 20; 
    melt_percentage(count) = pressure_pre;
    melt_diff(count) = pressure_diff;
    count = count+1;

end

    
%% plot figure 

figure('Position', [10 10 700 800]);
% plot P-T diagram
%line (melt_percentage, depth_z, 'Color', [0.8660 0.2740 0.1880],'LineWidth',2);
line (melt_diff, depth_z, 'Color', [0.8660 0.2740 0.1880],'LineWidth',2);
% line (melt_drop, depth_z, 'Color', [0.8660 0.2740 0.1880],'LineWidth',2);
% hold on; 
% line (melt, depth_z,'LineWidth',2, 'Color', [0 0.4470 0.7410]);
set(gca,'Ydir','reverse', 'FontSize',20);
%xlabel('Melting Volume (km^3/km)');
%xlabel('Difference in extracted melt (%)');
xlabel('Difference in melting (km^3/km)');
ylabel('Depth (km)')
%legend( {'lake level drop', 'no change of lake level'},'Location','southeast')

hold off;




%% export 
% Create a png

% %output_image = ['output/AfricaModels2022/Figure/',char(figname(end)),'_melting_history',num2str(nbegin), '-', num2str(nend)];
output_image = ['output/AfricaModels2022/Figure/extraction/',char(figname(end)),'_depth_meltV_plot_diff_',num2str(step), '_f', num2str(factor)];
print('-dpng','-r200',[output_image,'.png']);
 


end
