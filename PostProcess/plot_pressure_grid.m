function plot_pressure_grid(inpath,bpath, step1, step2, xlims,ylims,pr_type )
% This code will plot the lithostatic pressure rate or dynamic 
% pressrue rate change 

% Inpath = folder where grid_###.mat and markers_###.mat files  
% exist for first model 
% step1: step number for model 1
% step2: step number for model 2
% step2 should > step 1
% bpath, the path of bounds.mat
% xlims = [100 300];  % in km
% ylims = [-1 100]; % in km 
% gridsize : interpolated grid size i km 
% pr_type: the type of pressure, 1: lithostatic; 2 dynamic ; 3 total 

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
% e.g.: plot_pressure_grid('output/AfricaModels2022/Lake_model_extraction/Lake_drop_restart_71_dry',
% 'output/AfricaModels2022/Lake_model_extraction/reference_71_dry', 222,232, [0, 250], [20, 120], 1)


% default resolution
npx = 1500; %800;
npy = 750; %400;
tempc = 5;

cmap = 'jet(36)'; %colormap(parula(5))
sticky_layer = 20000;
yr2sec = 365.25*24*3600;
%time_diff: time span between these two step
time_diff = (step2-step1)*5; % the time step is 5 kyr



%% load grids
% Load files 1

infile = [inpath,'/grids_',num2str(step1),'.mat']; % as G
G1 = load(infile);
infile = [inpath,'/markers_',num2str(step1),'.mat']; % as G
load(infile);
%Load files 2
infile = [inpath,'/grids_',num2str(step2),'.mat']; % as G
G2 =load(infile); % grid of model 2

%% calculate lithostatic pressure and dynamic pressure 
% load bounds
infile = [bpath,'/bounds.mat']; % as G
B=load(infile);
%  modify boundary condition 
B.bright(:,4) = 1;
B.bright(:,1) = 0;
B.bleft(:,4) = 1;
B.bleft(:,1) = 0;  
B.btop(:, 3) =0;
B.bbottom(:, 3) = 0;

% calculate lithostatic pressure 
[x,y, pr1] = getlithostaticP(G1.G, B.bleft, B.bright, B.btop, B.bbottom, B.bintern);
[x,y, pr2] = getlithostaticP(G2.G, B.bleft, B.bright, B.btop, B.bbottom, B.bintern);

% lithostatic pressure rate in MPa
P_lith= 1e-6*(pr2- pr1)/time_diff; % MPa
% dynamic pressure rate in MPa
P_dyn= 1e-6*((G2.G.pr1- pr2)- (G1.G.pr1-pr1))/time_diff; % MPa
% total pressure in Mpa
P_tot = 1e-6*(G2.G.pr1 - G1.G.pr1)/time_diff; % MPa



%% interpolation 


%[Xq,Yq] = meshgrid(xlims(1):gridsize:xlims(2), ylims(1):gridsize:ylims(2));

%p1 = interp2(G.gridx(1:end-1)/1e3,G.gridy(1:end-1)/1e3, MGII1,Xq,Yq);
%p2 = interp2(G.gridx(1:end-1)/1e3,G.gridy(1:end-1)/1e3, MGII2,Xq,Yq);
%quiver (Xq, Yq, Vx,Vy,5, 'k');

%% Plot the material proporties difference
figure('Position', [10 10 1000 500]);


%time_diff = (step2 - step1)*5; % ky
if pr_type == 1
    p = pcolor(x/1e3, y/1e3, P_lith); % pressure change rate GPa/ky

    %hold on;
    %[con, h] = contourf(x/1e3, y/1e3, P_lith,20,'k','linewidth',1);  
    %[con, h] = contourf(x/1e3, y/1e3, P_lith,30);
    c = colorbar; 
    ylabel(c,'Lithostatic pressure change (MPa/kyr) ','FontSize',14);
    type= 'lithostatic';
elseif pr_type == 2
     p = pcolor(x/1e3, y/1e3, P_dyn); % pressure change rate GPa/ky
%     c = colorbar; 
%     ylabel(c,'Dynamic pressure change (GPa/kyr) ','FontSize',14);
%     hold on;
%     [con, h] = contour(x/1e3, y/1e3, P_dyn,10,'k','linewidth',1);
%    [con, h] = contourf(x/1e3, y/1e3, P_dyn,30);
    c = colorbar; 
    ylabel(c,'Dynamic pressure change (MPa/kyr) ','FontSize',14);
    type= 'dyanmic';
elseif pr_type == 3
     p = pcolor(x/1e3, y/1e3, P_tot); % pressure change rate GPa/ky
%     c = colorbar; 
%     ylabel(c,'Dynamic pressure change (GPa/kyr) ','FontSize',14);
%     hold on;
%     [con, h] = contour(x/1e3, y/1e3, P_dyn,10,'k','linewidth',1);
%    [con, h] = contourf(x/1e3, y/1e3, P_dyn,30);
    c = colorbar; 
    ylabel(c,'Total pressure change (MPa/kyr) ','FontSize',14);
    type= 'total';
else 
    ylabel('wrong input');
end
% setup the limits 

caxis([-5e-2, 5e-2]);
%caxis([-1e-1, 1e-1]);
colormap(cmap);

% plot lake 
mx = MX/1000;
% my = MY/1000 - sticky_layer;
my = MY/1000;

% Crop to desired size
k = mx >= xlims(1) & mx <= xlims(2) & my>= ylims(1) & my<= ylims(2);
mx1 = mx(k);
my1 = my(k);
mi = MI(k);
water = mi == 2;
trans_f = 0.3;
hold on;
wtr = [1, 1, 1];
marker_size = 3;
scatter(mx1(water),my1(water),marker_size,'MarkerFaceColor',wtr,'MarkerEdgeColor',wtr,...
    'MarkerFaceAlpha',trans_f,'MarkerEdgeAlpha',trans_f)

% plot melting triangle 
k2= MEXTC>0;
mx2 = mx(k2); my2 = my(k2);
% find the boundary 
k3 = boundary(mx2,my2);
plot(mx2(k3), my2(k3),'LineWidth', 2.5, Color = [1, 1, 1]); 
area = polyarea(mx2(k3),my2(k3));
disp(['area of melting triangle is ', area]);

ax=gca;
set(gca,'Ydir','reverse');
    % set up the scale for MGII
ax.FontSize = 14;
shading interp;


% title(['Step:',num2str(j),' Model time:',num2str(timesum*1e-6/(yr2sec)), ' Myr']);
xlabel('Horizon (km)','FontSize',14)
ylabel('Depth (km)', 'FontSize',14)
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
figname = split(inpath, '/');
%output_image = ['output/AfricaModels2022/Figure/extraction/',char(figname(end)),...
%     '_pressure_',type, '_rate_step',num2str(step1), '-',num2str(step2)];
%print('-dpng','-r200',[output_image,'.png'])


end

