function find_F(iframe)
% this is function plot hist and find F for ansthenosphere



load(['output/AfricaModels2022/Lake_model_extraction/reference_71/markers_',num2str(iframe),'.mat'])

%% plot figure 
figure();

pos = get(gcf,'Position');
width = pos(3) * 3;
height = pos(4) * 1.5;
set(gcf,'Position',[pos(1) pos(2) width height]);

i = MI == 10;
histogram(log10(META(i)), 30);
%histogram(log10(META(i)), 30, 'Binlimits', [19, 22]);
xlabel('Visocity (log)');
%set(gca, 'xscale','log');
set(gca,'FontSize',25);
%xlim ([]);
title(['Time = 5 Myr']);

%% export 
% % Create a png
output_image = ['output/AfricaModels2022/Figure/extraction/find_f_original']
print('-dpng','-r200',[output_image,'.png']);



end
