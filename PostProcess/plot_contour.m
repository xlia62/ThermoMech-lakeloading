function animate_plasticstrain(inpath,step,xlims,ylims, MPROP )
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

% xlims = [100 300];  % in km
% ylims = [-1 100]; % in km 
% npx = 1500;
% npy = 750;
npx = 1500;
npy = 750;
tempc = 5;

j = step;
cmap = jet; %colormap(parula(5))

% if nfreq == 1
%    k = [1 nfreq:nfreq:ntotal];
% else
%     k = [nfreq:nfreq:ntotal];
% end



yr2sec = 365.25*24*3600;

    
% Load files
infile = [inpath,'/markers_',num2str(j),'.mat'];
%load(infile,'MX','MY','MGII','MI', 'META');
load(infile);
infile = [inpath,'/grids_',num2str(j),'.mat']; % as G
load(infile);


 %% 

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
        unit = 'Cummulative Melt Fraction (%)';
    case 'MEXTC' % Cummulative Extracted Melt Fraction
        MGII = MEXTC;
        unit = 'Cummulative Extracted Melt Fraction';
    case 'MEXT' %  Extracted Melt Fraction
        MGII = MEXT;
        unit = 'Extracted Melt Fraction';
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

[mx,my,markstrain] = vis_markers((xlims-ox)*1000,(ylims-oy)*1000,npx,npy,mx,my,mgii);


mx = mx + ox*1000;
my = my + oy*1000;
    
if j == 1
    xlimits = [mx(1) mx(end)];
    ylimits = [my(1) my(end)];
end


%% Plot the material proporties
%clf
figure('Position', [10 10 1000 500]);

%pcolor(mx/1000,my/1000,-log10(abs(markstrain)))
plot(linspace(min(mx)/1e3, max(mx)/1e3, 100), 1e-3*water_lev*ones(100), ...
    'color', 'blue');
%legend({'lake level'}, 'Location','northwest');
hold on;

s = pcolor(mx/1000,my/1000,100*markstrain);
colormap(cmap);
c = colorbar; 
s.FaceColor = 'interp';

[con, h] = contour(G.gridx/1000,G.gridy/1000,G.tk1-273,[500 800 1000 1200 1400 1600],'white','linewidth',1);
clabel(con, h,  'Fontsize', 18, 'Color', 'white', 'labelspacing', 400);
ax=gca;
ax.FontSize = 18;
set(gca,'Ydir','reverse')
%title(['Step=',num2str(j),' Model time: ',num2str(timesum*1e-6/(yr2sec)), 'Myr']);
xlabel('Horizon (km)','FontSize',18)
ylabel('Depth (km)','FontSize',18)
ylabel(c,unit,'FontSize',18);
shading interp;
% set up the scale for MGII
% if MPROP == 'MGII'
%     caxis([0 10]);
%end
% if MPROP == 'MXM'
%     caxis([0.14 0.16]);
% end
%caxis([8.7 9.7]);

%legend( {'lake level'}, 'FontSize',16)
axis image
if ~isempty(xlims)
    xlim(xlims)
end
if ~isempty(ylims)
    ylim(ylims)
end
drawnow

    
% move the sticky air from figure
move = gca; 
tick = move.YTick - 20;
yticklabels({tick});

% Create a png
figname = split(inpath, '/');
output_image = ['output/AfricaModels2022/Figure/',char(figname(end)),'_step',num2str(step), '_2_', MPROP];
print('-dpng','-r200',[output_image,'.png'])

end

