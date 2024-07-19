function plot_particles(input_dir,iframe,xlims,ylims)

%input_dir = 'Africa2019/lake1200/Tmoho600_dT0';
output_image = [input_dir,'_frame',num2str(iframe)]

%nframes = 100;
sticky_layer = 20;
%xlims = [150 225];
%ylims = [-1 20];
marker_size = 10;


yr2sec= 365.25*24*3600;

% colors
c1 = [0.6350, 0.0780, 0.1840];
c2 = [0.8500, 0.3250, 0.0980];
c3 = [139 69 19]/255;
flt = [0.4940, 0.1840, 0.5560];
wtr = [0 0.447 0.7410];
sds = [0.9290 0.6940 0.1250];
mtl = [153 204 0]/255;
ast = [.47 .67 .19];


figure('Position',[66 395 1167 385],'color','white')

load([input_dir,'/markers_',num2str(iframe),'.mat'])
mx = MX/1000;
my = MY/1000 - sticky_layer;

k = mx >= xlims(1) & mx <= xlims(2) & my>= ylims(1) & my<= ylims(2);
mx = mx(k);
my = my(k);
mi = MI(k);

crust1 = mi == 4;
crust2 = mi == 6;
fault = mi == 11;
water = mi == 2;
seds = mi == 3;
crust3 = mi == 7;
mantle = mi == 9;
astheno = mi == 10;

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

axis image
ylim(ylims)
xlim(xlims)
set(gca,'Ydir','reverse')
set(gca,'FontSize',14)

title([num2str(timesum/yr2sec/1e6),' My'],'HorizontalAlignment','right','Position',[200 -1 0])

print('-dtiff','-r100',[output_image,'.tiff'])
