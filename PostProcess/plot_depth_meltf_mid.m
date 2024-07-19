function T=plot_depth_meltf_mid(inpath,step, maxdepth, nfreq)
% this is function plotting the melting fraction vs depth (%) for model 
% and reference model 

% input:
% Inpath = folder where grid_###.mat and markers_###.mat files exist
% step: the step in which melting volume should be calculated 
% nfreq: the numbers of depth interval of calculation 
% maxdepth for the plot in km 
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
%rho_solid = 3300;
%rho_melt = 2900;
%marker_per_grid = 36;
% factor = 1.01;
% output: 
% image of melting history(km^3/km)

% required:
% totalmelt_lake function

% main code:
%% load data and define indexes

depth = linspace(20, maxdepth, nfreq+1); % 20 km is the thickness of the sticky air 
interval =  (maxdepth-20)/nfreq; %depth interval
% Load files
infile = [inpath,'/markers_',num2str(step),'.mat'];
load(infile);

refername = replace(char(infile), 'drop', 'nochange'); 
% reference model 
ref = load(refername);
% As the width of model is 300, but it extends through time
% we calculate the mid line by
mid = 110e3; %0.5*(max(MX)- min(MX)); % m
% we need to set a tolerence for such mid zone
mid_tol = 5e3; %m 

% create a list for data
meltf = zeros(length(depth), 1);
meltf_r = zeros(length(depth), 1);
depth_z = zeros(length(depth), 1); % depth in km 
count = 1;

for d = depth
    
    depthlims = [d, d+interval]; % in km
    
    % now calculate each depth interval
    % find out the box in the depth interval with the mid zone 
    i = MY/1000 >= depthlims(1) & MY/1000 <= depthlims(2)& MX>= mid-mid_tol& MX <= mid+mid_tol;
    i_r = ref.MY/1000 >= depthlims(1) & ref.MY/1000 <= depthlims(2)& ref.MX>= mid-mid_tol& ref.MX <= mid+mid_tol;
    % cucmulative extracted melt fraction 
    mextc= mean(MEXTC(i));
    % for the refrence model 
    mextc_r = mean(ref.MEXTC(i_r));
    % write file
    meltf(count) = mextc*100; %%
    meltf_r(count) = mextc_r*100;%%
    depth_z(count) = d- 20; 
    count = count+1;


end

    
%% plot figure 

figure('Position', [10 10 700 800]);
% plot P-T diagram
%line (melt_percentage, depth_z, 'Color', [0.8660 0.2740 0.1880],'LineWidth',2);
line (meltf, depth_z, 'Color', [0.8660 0.2740 0.1880],'LineWidth',2);
hold on; 
line (meltf_r, depth_z,'LineWidth',2, 'Color', [0 0.4470 0.7410],'LineWidth', 2);
set(gca,'Ydir','reverse', 'FontSize',20);
%xlabel('Melting Volume (km^3/km)');
%xlabel('Difference in extracted melt (%)');
xlabel('Cumulative melt fraction (%)');
ylabel('Depth (km)')
legend( {'model with lake level change', 'no change of lake level'},'Location','southeast')

hold off;




%% export 
% Create a png
figname = split(inpath, '/'); 
% %output_image = ['output/AfricaModels2022/Figure/',char(figname(end)),'_melting_history',num2str(nbegin), '-', num2str(nend)];
output_image = ['output/AfricaModels2022/Figure/extraction/',char(figname(end)),'_depth_meltVf_plot_middle_',num2str(step)];
print('-dpng','-r200',[output_image,'.png']);
% %  


end
