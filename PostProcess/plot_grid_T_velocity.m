function plot_grid_T_velocity(inpath,step, xlims,ylims)
%% 

% Inpath = folder where grid_###.mat and markers_###.mat files exist


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
npx = 800;
npy = 400;
tempc = 5;

j = step;
cmap = 'jet'; %colormap(parula(5))

% if nfreq == 1
%    k = [1 nfreq:nfreq:ntotal];
% else
%     k = [nfreq:nfreq:ntotal];
% end



yr2sec = 365.25*24*3600;




%% 
% Load files 1
% infile = [inpath1,'/markers_',num2str(j),'.mat'];
% %load(infile,'MX','MY','MGII','MI', 'META');
% load(infile);
infile = [inpath,'/grids_',num2str(j),'.mat']; % as G
load(infile);


%% Plot the material proporties difference
%figure('Position', [10 10 1000 500]);
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

% p = pcolor(mx/1000,my/1000, markstrain1- markstrain2);
p = pcolor (G.gridx/1e3, G.gridy/1e3, G.tk1-273);
set(p,'facealpha',0.5);
%p = pcolor (G.gridx(2:end)/1e3, G.gridy(2:end)/1e3, MGII2);
%p = pcolor(mx/1000,my/1000, mgii);
% setup the limits 
%caxis([-1e8, 1e8]);
colormap(cmap);
c = colorbar; 
% ylabel(c,'Viscosity (10^{{\itx}} Pa s)','FontSize',20);
ylabel(c,['Temperature (', char(176),'C)'],'FontSize',16);

hold on
if ~isempty(tempc)
    [con, h] = contour(G.gridx/1000,G.gridy/1000,G.tk1-273,[500 800 1000 1200 1500],'white','linewidth',1);
    clabel(con, h,  'Fontsize', 14, 'Color', 'white', 'labelspacing', 400);
    %caxis([0 1450]);
end

%hold off
%colormap(cmap);

% plot velocity 
[Xq,Yq] = meshgrid(xlims(1):2:xlims(2), ylims(1):2:ylims(2));
% x =  400;
% y = 80;
%[Xq,Yq] = meshgrid(x,y);
[Vx, Vy] = velocityplot(Xq*1e3,Yq*1e3,G);
quiver (Xq, Yq, Vx,Vy,1, 'k');



shading interp;
%shading flat;

    
set(gca,'Ydir','reverse')
%title(['Step:',num2str(j),' Model time:',num2str(timesum*1e-6/(yr2sec)), ' Myr']);
xlabel('x, km')
ylabel('y, km')
%caxis([-1e7 1e7])
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
figname = split(inpath, '/');
output_image = ['output/AfricaModels2022/Figure/',char(figname(end)),'_velocity_T_',num2str(step)];
print('-dpng','-r200',[output_image,'.png'])


end

