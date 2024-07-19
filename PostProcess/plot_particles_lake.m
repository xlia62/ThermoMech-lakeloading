function output_image=plot_particles2(input_dir,iframe,xlims,ylims,marker_size)

% input_dir = 'Africa2019/lake1200/Tmoho600_dT0';
% iframe = 100;
% xlims = [150 225];  % in km
% ylims = [-1 20];    % in km
% sticky_layer = 20;  % in km

%marker_size = 20;
sticky_layer = 0; 


yr2sec= 365.25*24*3600;

% colors for markers
c1 = [0.6350, 0.0780, 0.1840];
c2 = [0.8500, 0.3250, 0.0980];
c3 = [139 69 19]/255;
flt = [0.4940, 0.1840, 0.5560];
wtr = [0 0.447 0.7410];
sds = [0.9290 0.6940 0.1250];
mtl = [153 204 0]/255;
ast = [.47 .67 .19];
% bslt = [210,105,30]/255;
bslt = [255,0,0]/255;
icc = [0 0.5 0.7410];

% Position and size of figure on computer screen
%figure('Position',[66 395 1167 385],'color','white')
figure('Position',[17 130 1262 489],'color','white')

% Load the marker file for the iframe
load([input_dir,'/markers_',num2str(iframe),'.mat'],'MX','MY','MI','timesum')
load([input_dir,'/grids_',num2str(iframe),'.mat'])

mx = MX/1000;
my = MY/1000 - sticky_layer;

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
basalt = mi == 12;  % Basalt
ice = mi == 18; %rift lake

scatter(mx(astheno),my(astheno),marker_size,'MarkerFaceColor',ast,'MarkerEdgeColor',ast,...
    'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)
hold on
scatter(mx(mantle),my(mantle),marker_size,'MarkerFaceColor',mtl,'MarkerEdgeColor',mtl,...
    'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)
scatter(mx(crust3),my(crust3),marker_size,'MarkerFaceColor',c3,'MarkerEdgeColor',c3,...
    'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);
scatter(mx(crust2),my(crust2),marker_size,'MarkerFaceColor',c2,'MarkerEdgeColor',c2,...
    'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)
scatter(mx(crust1),my(crust1),marker_size,'MarkerFaceColor',c1,'MarkerEdgeColor',c1,...
    'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)
scatter(mx(fault),my(fault),marker_size,'MarkerFaceColor',flt,'MarkerEdgeColor',flt,...
    'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)
scatter(mx(seds),my(seds),marker_size,'MarkerFaceColor',sds,'MarkerEdgeColor',sds,...
     'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)
scatter(mx(water),my(water),marker_size,'MarkerFaceColor',wtr,'MarkerEdgeColor',wtr,...
    'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)
% scatter(mx(basalt),my(basalt),marker_size,'MarkerFaceColor',bslt,'MarkerEdgeColor',bslt,...
%     'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)
% scatter(mx(ice),my(ice),marker_size,'MarkerFaceColor',icc,'MarkerEdgeColor',icc,...
%     'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)
% scatter(mx(mark),my(mark),marker_size,'MarkerFaceColor',icc,'MarkerEdgeColor',icc,...
%     'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)

% Dress the image
%axis image
ylim(ylims)
xlim(xlims)
xlabel('Horizontal position (km)')
ylabel('Depth (km)')
set(gca,'Ydir','reverse')
set(gca,'FontSize',20)


% move the sticky air from figure
move = gca; 
tick = move.YTick - 20;
yticklabels({tick});

%   Plot horizons
modeltime = timesum/yr2sec/1e6;
% nh = floor(modeltime)+1;
% for i = 1:nh
%     plot(gridt(1,:)/1000,H(i,:)/1000-sticky_layer,'linewidth',1,'color','b')
% end
plot(gridt(1,:)/1000,gridt(2,:)/1000-sticky_layer,'linewidth',1,'color','k')
% 
% title([num2str(modeltime),' My'],'HorizontalAlignment','right','Position',[200 -1 0])
% title(['Step:',num2str(iframe),' Model time',num2str(timesum*1e-6/(yr2sec)), ' Myr']);

%title([num2str(timesum/yr2sec/1e6),' My'],'HorizontalAlignment','right','Position',[200 -1 0])
title(['Step:',num2str(iframe),', Model time: ',num2str(timesum*1e-6/(yr2sec)), ' Myr']);

figname = split(input_dir, '/');
output_image = ['output/AfricaModels2022/Figure/extraction/',char(figname(end)),'_newframe_',num2str(iframe)];
print('-dpng','-r100',[output_image,'.png'])

% Create a tiff
% print('-dtiff','-r100',[output_image,'.tiff'])

end
% scatter(mx(water),my(water),marker_size,'MarkerFaceColor',wtr,'MarkerEdgeColor',wtr,...
%     'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)
