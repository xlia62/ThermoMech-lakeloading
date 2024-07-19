function plot_t_meltV(volume)
% this is function plot the melting history(km^3/km)

% input:
% volume: 1 volume, or 2: max melt fraction 
% melt_i choose to calculate melt or extract 1:MXM or 2:MEXTC



% output: 
% image of differental melting history
% required:
% totalmelt_lake function

% main code:
%% load data and define indexes

factor = 30; 

strong =  readtable('output/AfricaModels2022/Table/reference_weakcrust_melting_extract_history_1-71_dt_1_nomark_36_x_0-300_y_0-200.csv');
weak =  readtable('output/AfricaModels2022/Table/reference_strongcrust_melting_extract_history_1-71_dt_1_nomark_36_x_0-300_y_0-200.csv');
tp10 =  readtable('output/AfricaModels2022/Table/reference_Tp30_melting_extract_history_1-71_dt_1_nomark_36_x_0-300_y_0-200.csv');
lab65 =  readtable('output/AfricaModels2022/Table/reference_LAB65_melting_extract_history_1-71_dt_1_nomark_36_x_0-300_y_0-200.csv');
exten27 =  readtable('output/AfricaModels2022/Table/reference_v27_melting_extract_history_1-71_dt_1_nomark_36_x_0-300_y_0-200.csv');
ref = readtable('output/AfricaModels2022/Table/reference_71_melting_extract_history_1-71_dt_1_nomark_36_x_0-300_y_0-200.csv');
dry = readtable('output/AfricaModels2022/Table/reference_71_melting_extract_history_1-71_dt_1_nomark_36_x_0-300_y_0-200.csv');

figname = split('output/AfricaModels2022/Table/reference_dry_melting_extract_history_1-125_dt_1_nomark_36_x_0-300_y_0-200', '/'); 


name = char(figname(end));

%% plot figure 
figure();

pos = get(gcf,'Position');
width = pos(3) * 3;
height = pos(4) * 1.5;
set(gcf,'Position',[pos(1) pos(2) width height]);
%yyaxis left
if volume == 1 % plot melting volume 

elseif volume == 2  % plot melting percentage
    plot(strong.time_ky/1e3, factor*strong.max_mf, '-o','linewidth',4);
    hold on;
    plot(weak.time_ky/1e3, factor*weak.max_mf, '-o','linewidth',4);
    plot(tp10.time_ky/1e3, factor*tp10.max_mf, '-o','linewidth',4);
    plot(lab65.time_ky/1e3, factor*lab65.max_mf, '-o','linewidth',4);
    plot(exten27.time_ky/1e3, factor*exten27.max_mf, '-o','linewidth',4);
    plot(ref.time_ky/1e3, factor*ref.max_mf, '-o','linewidth',4);
    plot(dry.time_ky/1e3, factor*dry.max_mf, '-o','linewidth',4);
    ylabel({'Melt fraction (%)'});

else
    disp('please choose plot melting volume or percentage')
end 

grid on;
xlabel('Model time (Myr)');



% yyaxis right

% plot (ref.time_ky, ref.water_mass_kg,'linewidth',4, 'color', [0 0.4470 0.7410]);
% % plot ((T_600.time_ky -T_600.time_ky(1))/1e3, T_600.water_mass_kg,'linewidth',2);
% ylabel('Water mass (kg/km)');

legend({'model with strong crust','model with weak crust', 'model with Tp = 1353°C','model with TL = 105 km', ...
    'half extension 3 mm/yr', 'reference model', 'dry mantle'}, 'Location','northwest');

set(gca,'FontSize',25);
% ax = gca;
% ax.YAxis(1).Color = [0.8500 0.3250 0.0980];
% ax.YAxis(2).Color = [0 0.4470 0.7410];

hold off;


%% export 
% % Create a png
% output_image = ['output/AfricaModels2022/Figure/extraction/',name(1:end-4),...
%  '_melt_fraction_reference'];
% print('-dpng','-r200',[output_image,'.png']);

% time_ky = time*1e3;
% melt_v = melt_drop'; %km^2/km
% melt_refer = melt; %km^2/km
% water_mass_kg = water; %kg
% max_mf = melt_fraction;  % max accumulated melt fraction
% T = table(time_ky, melt_v, melt_refer, water_mass_kg, max_mf);
% output_table = ['output/AfricaModels2022/Table/',char(figname(end)),'_',...
%     num2str(type),'_history_',num2str(nbegin), '-', num2str(nend), '_dt_', num2str(nfreq)];
% writetable(T, [output_table, '.csv']); 


end
