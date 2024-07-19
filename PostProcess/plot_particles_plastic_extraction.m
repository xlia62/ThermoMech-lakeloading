function output_image=plot_particles_plastic(input_dir,iframe,threshold,thresmax, type)
% type 1 : all the basin 
% type 2 : border fault 
% type 3 : secondary border fault 



% input_dir = 'Africa2019/lake1200/Tmoho600_dT0';
% iframe = 100;
% xlims = [150 225];  % in km
% ylims = [-1 20];    % in km
% threshold of the plastic strain
sticky_layer = 20;  % in km
marker_size = 5;
% npx = 1000;
% npy = 500;

cmap = jet;

yr2sec= 365.25*24*3600;

% colors for markers
gray = [0, 0, 0];
c1 = gray +0.05*12 ;
c2 =  gray +0.05*17 ;
c3 =  gray +0.05*14 ;
%flt =  gray +0.05*4 ;
wtr =  gray +0.05*8 ;
mtl =  gray +0.05*8 ;
ast =  gray +0.05*5 ;


% Position and size of figure on computer screen
%figure('Position',[66 395 1167 385],'color','white')
%figure('Position',[17 130 1262 489])
%figure('Position',[17 130 1362 489])
figure();
pos = get(gcf,'Position');

if type == 1
    xlims = [50 200];
    ylims = [15 30];
    width = pos(3) * 3;
    height = pos(4) * 1.5;
    npx = 3000;
    npy = 300;
elseif type == 2
    xlims = [110 160];
    ylims = [15 30];  
    width = pos(3) * 2.5;
    height = pos(4) * 1.5;
    npx = 1000;
    npy = 300;
elseif type == 3
    xlims = [55 110]; % 55
    ylims = [15 40];  % 25
    width = pos(3) * 2.5;
    height = pos(4) * 1.5;
    npx = 1100;
    npy = 500;
else
    disp('wrong input');
end


set(gcf,'Position',[pos(1) pos(2) width height]);


% Load the marker file for the iframe
load([input_dir,'/markers_',num2str(iframe),'.mat'])
load([input_dir,'/grids_',num2str(iframe),'.mat'])

mx = MX/1000;
% my = MY/1000 - sticky_layer;
my = MY/1000;

% Crop to desired size
k = mx >= xlims(1) & mx <= xlims(2) & my>= ylims(1) & my<= ylims(2);
mx = mx(k);
my = my(k);
mi = MI(k);

% Mapping of markers to rocktypes or fluids
crust1 = mi == 4;  % upper 1
crust2 = mi == 6;  % upper 2
mark = mi == 5;  % marker of upper 
fault = mi == 11;  % Weak-zone
water = mi == 2;
seds = mi == 3;
crust3 = mi == 7;  % lower crust
mantle = mi == 9;  % Mantle lithosphere
astheno = mi == 10; 
ice = mi == 18; %rift lake
trans_f = 0.3;

scatter(mx(astheno),my(astheno),marker_size,'MarkerFaceColor',ast,'MarkerEdgeColor',ast,...
    'MarkerFaceAlpha',trans_f,'MarkerEdgeAlpha',trans_f)
hold on
scatter(mx(mantle),my(mantle),marker_size,'MarkerFaceColor',mtl,'MarkerEdgeColor',mtl,...
    'MarkerFaceAlpha',trans_f,'MarkerEdgeAlpha',trans_f)
scatter(mx(crust3),my(crust3),marker_size,'MarkerFaceColor',c3,'MarkerEdgeColor',c3,...
    'MarkerFaceAlpha',trans_f,'MarkerEdgeAlpha',trans_f);
scatter(mx(crust2),my(crust2),marker_size,'MarkerFaceColor',c2,'MarkerEdgeColor',c2,...
    'MarkerFaceAlpha',trans_f,'MarkerEdgeAlpha',trans_f)
scatter(mx(crust1),my(crust1),marker_size,'MarkerFaceColor',c1,'MarkerEdgeColor',c1,...
    'MarkerFaceAlpha',trans_f,'MarkerEdgeAlpha',trans_f)
% scatter(mx(fault),my(fault),marker_size,'MarkerFaceColor',flt,'MarkerEdgeColor',flt,...
%     'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)
% scatter(mx(seds),my(seds),marker_size,'MarkerFaceColor',sds,'MarkerEdgeColor',sds,...
%     'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2, )
scatter(mx(water),my(water),marker_size,'MarkerFaceColor',wtr,'MarkerEdgeColor',wtr,...
    'MarkerFaceAlpha',trans_f,'MarkerEdgeAlpha',trans_f)
% scatter(mx(ice),my(ice),marker_size,'MarkerFaceColor',icc,'MarkerEdgeColor',icc,...
%     'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2, )
% scatter(mx(mark),my(mark),marker_size,'MarkerFaceColor',icc,'MarkerEdgeColor',icc,...
%     'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2, )


% plot strain 
load([input_dir,'/markers_',num2str(iframe),'.mat'])
load([input_dir,'/grids_',num2str(iframe),'.mat'])
MGII = MGII;
unit = {'Accumulated plastic strain (lg)'};
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
    
s = pcolor(mx/1000,my/1000,markstrain);
% s = pcolor(mx/1000,my/1000,log10(abs(markstrain)));
%s = pcolor(mx/1000,my/1000,log10(markstrain));
colormap(cmap);
c = colorbar; 
s.FaceColor = 'interp';
    
% Dress the image
%axis image
ylim(ylims)
xlim(xlims)
xlabel('Horizon (km)')
ylabel('Depth (km)')
set(gca,'Ydir','reverse')
set(s, 'EdgeColor', 'none');
set(gca,'FontSize',25)
%set(s, 'AlphaData', markstrain > threshold, 'FaceAlpha',  'interp')  %'FaceAlpha',.3,'EdgeAlpha',.3
%set(s, 'AlphaData', markstrain > log10(threshold), 'FaceAlpha',  'interp') 
%set(s, 'FaceAlpha',.8,'EdgeAlpha',.8)
ylabel(c,[unit],'FontSize',25);
%caxis([log10(threshold) log10(thresmax)])
caxis([threshold thresmax])

%   Plot horizons
modeltime = timesum/yr2sec/1e6;
% move the sticky air from figure
move = gca; 
tick = move.YTick - 20;
yticklabels({tick});

%title([num2str(timesum/yr2sec/1e6),' My'],'HorizontalAlignment','right','Position',[200 -1 0])
title(['Step:',num2str(iframe),', Model time: ',num2str(timesum*1e-6/(yr2sec)), ' Myr']);

figname = split(input_dir, '/');
output_image = ['output/AfricaModels2022/Figure/extraction/',char(figname(end)),'_frame',num2str(iframe),...
    '_P_log_', num2str(threshold), '-', num2str(thresmax),'_type', num2str(type), '_'];
print('-dpng','-r100',[output_image,'.png'])

% Create a tiff
% print('-dtiff','-r100',[output_image,'.tiff'])

end
% scatter(mx(water),my(water),marker_size,'MarkerFaceColor',wtr,'MarkerEdgeColor',wtr,...
%     'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)