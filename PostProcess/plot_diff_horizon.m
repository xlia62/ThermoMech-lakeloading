function output_image=plot_particles2(input_dir1,input_dir2, iframe,xlims,ylims)
% this function plot the difference between two horizons 

% input:
% input_dir12 = folder where grid_###.mat and markers_###.mat files exist
% iframe = step of the model, note that frame should be the same for
% comparision
% xlims,ylims = constrain caclulations to x-limits and y-limits of model, in km
% marker size: the size of marker 

% output: 
% image of different horizons 

% varibles
marker_size = 15;
sticky_layer = 20; % km, thickness of sticky layer
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
bslt = [47 79 79]/255;
icc = [0 0.5 0.7410];

% main code:
%% load data
% Load the marker file1 for the iframe
load([input_dir1,'/markers_',num2str(iframe),'.mat'],'MX','MY','MI','timesum')
load([input_dir1,'/grids_',num2str(iframe),'.mat'])

mx1 = MX/1000;
my1 = MY/1000 - sticky_layer;

% Crop to desired size
k = mx1 >= xlims(1) & mx1 <= xlims(2) & my1>= ylims(1) & my1<= ylims(2);
mx1 = mx1(k);
my1 = my1(k);
mi1 = MI(k);
gridt1 = gridt;
% Mapping of markers to rocktypes or fluids
marker1 = mi1 == 5;  % marker betwen upper and lower 
water = mi1 == 2;

% Load the marker file2 for the iframe
load([input_dir2,'/markers_',num2str(iframe),'.mat'],'MX','MY','MI','timesum')
load([input_dir2,'/grids_',num2str(iframe),'.mat'])

mx2 = MX/1000;
my2 = MY/1000 - sticky_layer;

% Crop to desired size
k = mx2 >= xlims(1) & mx2 <= xlims(2) & my2>= ylims(1) & my2<= ylims(2);
mx2 = mx2(k);
my2 = my2(k);
mi2 = MI(k);
gridt2 = gridt;
% Mapping of markers to rocktypes or fluids
marker2 = mi2 == 5;  % marker betwen upper and lower 

% calculate the topography elevation difference
gridt_diff = gridt1- gridt2;

%% plot figure 
figure('Position',[17 130 1262 489],'color','white')
%figure;
clf;
% scatter(mx1(marker1),my1(marker1),marker_size,'MarkerFaceColor',c2,'MarkerEdgeColor',c2,...
%     'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)
% hold on;
% scatter(mx2(marker2),my2(marker2),marker_size,'MarkerFaceColor',bslt,'MarkerEdgeColor',bslt,...
%     'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)
yyaxis left
scatter(mx2(marker2),my1(marker1)-my2(marker2),marker_size,'MarkerFaceColor',[0.8500 0.3250 0.0980],'MarkerEdgeColor',[0.8500 0.3250 0.0980],...
    'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)
hold on;
scatter(gridt1(1,:)/1000, gridt_diff(2,:)/1000,marker_size,'MarkerFaceColor',[0.4660 0.6740 0.1880],'MarkerEdgeColor',[0.4660 0.6740 0.1880],...
    'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)

% Dress the image
% axis image
% %ylim(ylims)
xlim(xlims)
xlabel('X (km)')
ylabel('Elevation difference (km)')
set(gca,'FontSize',16)

%grid on;

yyaxis right

% scatter(mx1(water), -my1(water)+300,0.5,'MarkerFaceColor',icc,'MarkerEdgeColor',icc,...
%     'MarkerFaceAlpha',.01,'MarkerEdgeAlpha',.01)
% plot (gridt2(1,:)/1000, (-gridt2(2,:)+ 320e3)/1000, 'k','LineStyle','-','linewidth',1);
scatter(mx1(water), -my1(water)+300,5,'MarkerFaceColor',icc,'MarkerEdgeColor',icc,...
    'MarkerFaceAlpha',.1,'MarkerEdgeAlpha',.1)
plot (gridt2(1,:)/1000, (-gridt2(2,:)+ 320e3)/1000, 'k','LineStyle','-','linewidth',1);

xlim(xlims)
ylim([200 310])
xlabel('X (km)')
ylabel('topography (km)')
legend('elevation difference of upper/lower crust', 'elevation difference of topography', 'lake','topography',...
    'Location','southeast');

modeltime = timesum/yr2sec/1e6;
title(['Step:',num2str(iframe),', Model time: ',num2str(timesum*1e-6/(yr2sec)), ' Myr']);

ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';

figname = split(input_dir1, '/');
output_image = ['output/AfricaModels2022/Figure/',char(figname(end)),'_horizon_diff_frame',num2str(iframe)];
print('-dpng','-r100',[output_image,'.png'])

% Create a tiff
% print('-dtiff','-r100',[output_image,'.tiff'])

end
% scatter(mx(water),my(water),marker_size,'MarkerFaceColor',wtr,'MarkerEdgeColor',wtr,...
%     'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)