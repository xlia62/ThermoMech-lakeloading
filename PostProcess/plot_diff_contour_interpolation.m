function plot_diff_contour(inpath1,inpath2,step1,step2, gridsize, xlims,ylims, MPROP )
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
% default resolution 
npx = 100;
npy = 50;
tempc = 5;

% gridsize: interpolated grid size in km (1km )

cmap = 'jet'; %colormap(parula(5))

% if nfreq == 1
%    k = [1 nfreq:nfreq:ntotal];
% else
%     k = [nfreq:nfreq:ntotal];
% end



yr2sec = 365.25*24*3600;




%% 
% Load files 1

j = step1;
infile = [inpath1,'/markers_',num2str(j),'.mat'];
%load(infile,'MX','MY','MGII','MI', 'META');
load(infile);
infile = [inpath1,'/grids_',num2str(j),'.mat']; % as G
load(infile);


switch MPROP
    case 'META' % viscosity, Pa s
        MGII = META;
        unit = 'viscosity (Pa s)';
    case 'MGII' % Accumulated plastic strain
        MGII = MGII;
        unit = 'Accumulated plastic strain';
    case 'MBII' % Accumulated bulk strain
        MGII = MBII;
        unit = 'Accumulated bulk strain';
    case 'MPR' % pressure
        MGII = MPR;
        unit = 'Pressure (Pa)';
    case 'MTK' % Temperature, K
        MGII = MTK;
        unit = 'Temperature (K)';
    case 'MXM' % Cummulative Melt Fraction   
        MGII = MXM;
        unit = 'Cummulative Melt Fraction';
    case 'MEXTC' % Cummulative Extracted Melt Fraction
        MGII = MEXTC;
        unit = 'Cummulative Extracted Melt Fraction';
    case 'MSXY' % Cummulative Extracted Melt Fraction
        MGII = MSXY;
        unit = 'shear stress (Pa)';
    case 'MSXX' % Cummulative Extracted Melt Fraction
        MGII = MSXX;
        unit = 'normal stress (Pa)';
end    


% Replace Air or Water with NaN
i = MI == 1 | MI == 2;
MGII(i) = NaN;
    
% select based on x and y limits
if ~isempty(xlims)
    i = MX/1000 >= xlims(1) & MX/1000 <= xlims(2);
    ox = xlims(1);
    mx = MX(i) - ox*1000;
    my = MY(i);
    mgii = MGII(i);
else
    xlims = [G.gridx(1) G.gridx(end)]/1000;
    ox = 0;
    mx = MX;
    my = MY;
    mgii = MGII;
end
if ~isempty(ylims)
    i = my/1000 >= ylims(1) & my/1000 <= ylims(2);
    oy = ylims(1);
    my = my(i) - oy*1000;
    mx = mx(i);
    mgii = mgii(i);
else
    ylims = [G.gridy(1) G.gridy(end)]/1000;
    oy = 0;
    my = my;
end


%mgii = max(mgii,1e-5);
% 
[mx,my,markstrain1] = vis_markers((xlims-ox)*1000,(ylims-oy)*1000,npx,npy,mx,my,mgii);

mx = mx + ox*1000;
my = my + oy*1000;
    
if j == 1
    xlimits = [mx(1) mx(end)];
    ylimits = [my(1) my(end)];
end
    

%% Load files 2

j2 = step2;
infile = [inpath2,'/markers_',num2str(j2),'.mat'];
%load(infile,'MX','MY','MGII','MI', 'META');
load(infile);
infile = [inpath2,'/grids_',num2str(j2),'.mat']; % as G
load(infile);

switch MPROP
    case 'META' % viscosity, Pa s
        MGII = META;
        unit = 'viscosity (Pa s)';
    case 'MGII' % Accumulated plastic strain
        MGII = MGII;
        unit = 'Accumulated plastic strain';
    case 'MBII' % Accumulated bulk strain
        MGII = MBII;
        unit = 'Accumulated bulk strain';
    case 'MPR' % pressure
        MGII = MPR;
        unit = 'Pressure (Pa/yr)';
    case 'MTK' % Temperature, K
        MGII = MTK;
        unit = 'Temperature (K)';
    case 'MXM' % Cummulative Melt Fraction   
        MGII = MXM;
        unit = 'Cummulative Melt Fraction';
    case 'MEXTC' % Cummulative Extracted Melt Fraction
        MGII = MEXTC;
        unit = 'Cummulative Extracted Melt Fraction';
    case 'MSXY' % Cummulative Extracted Melt Fraction
        MGII = MSXY;
        unit = 'shear stress (Pa)';
    case 'MSXX' % Cummulative Extracted Melt Fraction
        MGII = MSXX;
        unit = 'normal stress (Pa)';
end    


% Replace Air or Water with NaN
i = MI == 1 | MI == 2;
MGII(i) = NaN;
    
% select based on x and y limits
if ~isempty(xlims)
    i = MX/1000 >= xlims(1) & MX/1000 <= xlims(2);
    ox = xlims(1);
    mx = MX(i) - ox*1000;
    my = MY(i);
    mgii = MGII(i);
else
    xlims = [G.gridx(1) G.gridx(end)]/1000;
    ox = 0;
    mx = MX;
    my = MY;
    mgii = MGII;
end
if ~isempty(ylims)
    i = my/1000 >= ylims(1) & my/1000 <= ylims(2);
    oy = ylims(1);
    my = my(i) - oy*1000;
    mx = mx(i);
    mgii = mgii(i);
else
    ylims = [G.gridy(1) G.gridy(end)]/1000;
    oy = 0;
    my = my;
end
    


[mx,my,markstrain2] = vis_markers((xlims-ox)*1000,(ylims-oy)*1000,npx,npy,mx,my,mgii);

mx = mx + ox*1000;
my = my + oy*1000;
    
if j2 == 1
    xlimits = [mx(1) mx(end)];
    ylimits = [my(1) my(end)];
end
%% interploation 
% [Xq,Yq] = meshgrid(xlims(1):gridsize:xlims(2), ylims(1):gridsize:ylims(2));
% 
% p1 = interp2(G.gridx(1:end-1)/1e3,G.gridy(1:end-1)/1e3, markstrain1',Xq,Yq);
% p2 = interp2(G.gridx(1:end-1)/1e3,G.gridy(1:end-1)/1e3, markstrain2',Xq,Yq);
%% Plot the material proporties difference
%clf
figure('Position', [10 10 1000 500]);
%pcolor(mx/1000,my/1000,-log10(abs(markstrain)))
colors1 = [.7 .7 .7; .4 .4 .4; .005 .005 .005];
n = 100;
m = size(colors1,1);
t0 = linspace(0,1,m)';
t = linspace(0,1,n)';
r = interp1(t0,colors1(:,1),t);
g = interp1(t0,colors1(:,2),t);
b = interp1(t0,colors1(:,3),t);
cmap2 = [r,g,b];
cmap3 = [0 1 0; 
         1 0 0;
         0 0 1];

p = pcolor(mx/1000,my/1000, (markstrain1- markstrain2)/5000);% the interval is 5kyr so the unit is Pa/yr
% setup the limits 
%caxis([-0.05, 0.05]);
colormap(cmap);
c = colorbar; 
% ylabel(c,'Viscosity (10^{{\itx}} Pa s)','FontSize',20);
ylabel(c,unit,'FontSize',16);

hold on
% if ~isempty(tempc)
%     h=contour(G.gridx/1000,G.gridy/1000,G.tk1-273,tempc,'k','linewidth',1);
% end
% hold off
% colormap(cmap);




shading interp;
%shading flat;

    
set(gca,'Ydir','reverse')
%title(['Step=',num2str(j),' Myr=',num2str(timesum*1e-6/(yr2sec))]);
xlabel('x, km')
ylabel('y, km')
set(gca,'FontSize',18)
%set(gca,'clim',[-2 1])

%end
caxis([-2e2 2e2]);
set(gca,'FontSize',18);
xlabel('Horizontal distance (km)')
ylabel('Depth (km)')
axis image
if ~isempty(xlims)
    xlim(xlims)
end
if ~isempty(ylims)
    ylim(ylims)
end
%drawnow
move = gca; 
tick = move.YTick - 20;
yticklabels({tick});
    
% if ~isempty(tempc)
%     clabel(C,h,'manual','rotation',0)
% end



% Create a png
% figname = split(inpath1, '/');
% output_image = ['output/AfricaModels2022/Figure/extraction/',char(figname(end)),...
%     '_diff_contour_interpolation_step',num2str(step1), '-',num2str(step2),'_', MPROP];
% print('-dpng','-r200',[output_image,'.png'])


%% histogram 
% left = -1e-2;
% right = 1e-2;
% edges = linspace(left, right, 50);
% figure; 
% histogram(markstrain1- markstrain2, 'BinEdges',edges)
% % set(gca, 'YScale', 'log')
% grid on;
% xlim([left, right]);
% xlabel(unit, 'FontSize', 14);
% ylabel('Bin Count', 'FontSize', 14);
% title('Histogram of Data', 'FontSize', 14);
% output_image = ['output/AfricaModels2022/Figure/',char(figname(end)),'_diff_contour_histogram_step',num2str(step), '_', MPROP];
% print('-dpng','-r200',[output_image,'.png'])

end
