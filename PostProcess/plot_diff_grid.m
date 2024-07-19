function animate_plasticstrain(inpath1,inpath2,step1, step2, gridsize, xlims,ylims, MPROP )
%% 

% Inpath = folder where grid_###.mat and markers_###.mat files exist

% nfreq  = frequency of time stpes when files are saved
% ntotal = total number of time steps
% npx,npy = number of image pixels in x and y direction
% cmap = colormap
% tempc = temperature contours, ommit adding contours if empty
% inpath = 'output/AfricaModels2022/fastscape_line2';
% xlims = [0 400];  % in km
% ylims = [-1 200]; % in km 
% MPROP = META;
%step =200;
% inpath = 'utput/AfricaModels2022/lake_decrease_faultzone'
% step = 55
% xlims = [100 300];  % in km
% ylims = [-1 100]; % in km 
% gridsize : interpolated grid size i km 

% G.etas1 = zeros(ynum,xnum);       % Viscosity for shear stress
% G.etan1 = zeros(ynum-1,xnum-1);   % Viscosity for normal stress
% G.mus1 = zeros(ynum,xnum);        % Shear modulus for shear stress
% G.mun1 = zeros(ynum-1,xnum-1);    % Shear modulus for normal stress
% G.sxy1 = zeros(ynum,xnum);        % Shear stress
% G.sxx1 = zeros(ynum-1,xnum-1);    % Normal stress
% G.rho1 = zeros(ynum,xnum);        % Density
% G.tk1 = zeros(ynum,xnum);         % Old temperature
% G.tk2 = zeros(ynum,xnum);                        % New temperature
% G.rhocp1 = zeros(ynum,xnum);      % RHO*Cp (for temperature equation)
% G.kt1 = zeros(ynum,xnum);         % Thermal conductivity
% G.hr1 = zeros(ynum,xnum);         % Radiogenic heat production
% G.ha1 = zeros(ynum,xnum);         % Adiabatic heat production/consuming
% G.hph = zeros(ynum,xnum);         % Latent heat from phase changes
% G.he = zeros(ynum,xnum);          % Adiabat/Shear heating efficiency

% default resolution
npx = 1500; %800;
npy = 750; %400;
tempc = 5;
j = step1;
h = step2;
cmap = 'jet(4096)'; %colormap(parula(5))
sticky_layer = 20000;
yr2sec = 365.25*24*3600;




%% 
% Load files 1

infile = [inpath1,'/grids_',num2str(j),'.mat']; % as G
load(infile);

switch MPROP
    case 'sxy1' 
        MGII1 = G.sxy1;
        unit = 'Shear stress (Pa)';
    case 'sxx1'
        MGII1 = G.sxx1;
        unit = 'Normal stress (Pa)';
    case 'rho1' % Accumulated bulk strain
        MGII1 = G.rho1;
        unit = 'Density (km/m^3)';
    case 'pr1' % pressure
        MGII1 = G.pr1;
        unit = 'Pressure (Pa)';
    case 'tk1' % Temperature, K
        MGII1 = G.tk1;
        unit = 'Temperature (K)';
    case 'eii' % bulk strain   
        MGII1 = G.eii;
        unit = 'Grid strain rates';
end    

  
%% Load files 2
infile = [inpath2,'/grids_',num2str(h),'.mat']; % as G
load(infile);


switch MPROP
    case 'sxy1' 
        MGII2 = G.sxy1;
        unit = 'Shear stress (Pa)';
    case 'sxx1'
        MGII2 = G.sxx1;
        unit = 'Normal stress (Pa)';
    case 'rho1' % Accumulated bulk strain
        MGII2 = G.rho1;
        unit = 'Density (km/m^3)';
    case 'pr1' % pressure
        MGII2 = G.pr1;
        unit = 'Pressure (Pa)';
    case 'tk1' % Temperature, K
        MGII2 = G.tk1;
        unit = 'Temperature (K)';
    case 'eii' % bulk strain   
        MGII2 = G.eii;
        unit = 'Grid strain rates';
end    

%% interpolation 


[Xq,Yq] = meshgrid(xlims(1):gridsize:xlims(2), ylims(1):gridsize:ylims(2));

p1 = interp2(G.gridx(1:end-1)/1e3,G.gridy(1:end-1)/1e3, MGII1,Xq,Yq);
p2 = interp2(G.gridx(1:end-1)/1e3,G.gridy(1:end-1)/1e3, MGII2,Xq,Yq);
%quiver (Xq, Yq, Vx,Vy,5, 'k');

%% Plot the material proporties difference
figure('Position', [10 10 1000 500]);

%pcolor(mx/1000,my/1000,-log10(abs(markstrain)))
% colors1 = [.7 .7 .7; .4 .4 .4; .005 .005 .005];
% n = 100;
% m = size(colors1,1);
% t0 = linspace(0,1,m)';
% t = linspace(0,1,n)';
% r = interp1(t0,colors1(:,1),t);
% g = interp1(t0,colors1(:,2),t);
% b = interp1(t0,colors1(:,3),t);
% cmap2 = [r,g,b];
% cmap3 = [0 1 0; 
%          1 0 0;
%          0 0 1];
% output time step = 20 ky
time_diff = (step2 - step1)*5; % ky
p = pcolor(Xq, Yq, (p2- p1)/time_diff); % pressure change rate P/ky
% p = pcolor (G.gridx(2:end)/1e3, G.gridy(2:end)/1e3, MGII1-MGII2);

%p = pcolor (G.gridx(2:end)/1e3, G.gridy(2:end)/1e3, MGII2);
%p = pcolor(mx/1000,my/1000, mgii);
% setup the limits 
%caxis([-1e8, 1e8]);
colormap(cmap);
c = colorbar; 
% ylabel(c,'Viscosity (10^{{\itx}} Pa s)','FontSize',20);
ylabel(c,'rate of pressure change (Pa/ky) ','FontSize',13);


hold on
if ~isempty(tempc)
    [con, h] = contour(G.gridx/1000,G.gridy/1000,G.tk1-273,[500 800 1000 1200 1500],'k','linewidth',1);
    clabel(con, h,  'Fontsize', 14, 'Color', 'k', 'labelspacing', 400);
    %caxis([0 1450]);
end

shading interp;
%shading flat;

    
set(gca,'Ydir','reverse')
% title(['Step:',num2str(j),' Model time:',num2str(timesum*1e-6/(yr2sec)), ' Myr']);
xlabel('Horizon (km)')
ylabel('Depth (km)')
%caxis([-5e5 5e5])
%set(gca,'clim',[-1e-13 1e-13])

%end


axis image
if ~isempty(xlims)
    xlim(xlims)
end
if ~isempty(ylims)
    ylim(ylims)
end
%drawnow
    
% move the sticky air from figure
move = gca; 
tick = move.YTick - 20;
yticklabels({tick});

% Create a png
figname = split(inpath1, '/');
output_image = ['output/AfricaModels2022/Figure/',char(figname(end)),...
    '_diff_grid_step',num2str(step1), '-',num2str(step2),'_gridsize_',num2str(gridsize),'_', MPROP];
print('-dpng','-r200',[output_image,'.png'])


end

