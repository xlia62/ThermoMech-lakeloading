function plot_t_meltV(volume)
% this is function plot the melting history(km^3/km)

% input:
% volume: 1 volume, or 2: percentage 
% melt_i choose to calculate melt or extract 1:MXM or 2:MEXTC



% output: 
% image of differental melting history
% required:
% totalmelt_lake function

% main code:
%% load data and define indexes
factor = -1;  % this is for the case where strengh profile is wrong. 
T_1500 =  readtable('output/AfricaModels2022/Table/Lake_drop_restart_71_1500m_dry_melting_extract_history_197-253_dt_1_nomark_36_x_0-300_y_0-200.csv');
T_1000 =  readtable('output/AfricaModels2022/Table/Lake_drop_restart_71_dry_melting_extract_history_197-253_dt_1_nomark_36_x_0-300_y_0-200.csv');
T_600 =  readtable('output/AfricaModels2022/Table/Lake_drop_restart_71_600m_dry_melting_extract_history_197-253_dt_1_nomark_36_x_0-300_y_0-200.csv');
T_300 =  readtable('output/AfricaModels2022/Table/Lake_drop_restart_71_300m_dry_melting_extract_history_197-253_dt_1_nomark_36_x_0-300_y_0-200.csv');
T_600 =  readtable('output/AfricaModels2022/Table/Lake_drop_restart_71_600m_melting_extract_history_367-435_dt_1_nomark_36_x_0-300_y_0-200.csv');
T_300 =  readtable('output/AfricaModels2022/Table/Lake_drop_restart_71_300m_melting_extract_history_367-435_dt_1_nomark_36_x_0-300_y_0-200.csv');

figname = split('output/AfricaModels2022/Table/Lake_drop_restart_71_300m_dry_melting_extract_history_197-253_dt_1_nomark_36_x_0-300_y_0-200.csv', '/'); 

T_1000.percentage = (T_1000.melt_v - T_1000.melt_refer)./(T_1000.melt_refer) ; 
T_600.percentage = (T_600.melt_v - T_600.melt_refer)./(T_600.melt_refer) ; 
T_300.percentage = (T_300.melt_v - T_300.melt_refer)./(T_300.melt_refer) ; 

% consider if ploting every 2 element 
% note that the total number is odd
melt_1500 = mean(reshape(T_1500.melt_v(2:end), 2, []));
melt_300 = mean(reshape(T_300.melt_v(2:end), 2, []));
melt_600 = mean(reshape(T_600.melt_v(2:end), 2, []));
melt_1000 = mean(reshape(T_1000.melt_v(2:end), 2, []));
melt_refer2 = mean(reshape(T_1500.melt_refer(2:end), 2, []));
melt_refer1 = mean(reshape(T_300.melt_refer(2:end), 2, []));

time1 = mean(reshape((T_300.time_ky(2:end) - T_300.time_ky(1)), 2, []));
time2 = mean(reshape((T_1500.time_ky(2:end) - T_1500.time_ky(1)), 2, []));
len = length(time2);
trend = linspace(0, 4.5, len);
percentage_1500 = 100*(melt_1500 - melt_refer2)./melt_refer2 +trend;
percentage_300 = 100*(melt_300 - melt_refer1)./melt_refer1;
percentage_600 = 100*(melt_600 - melt_refer1)./melt_refer1;
percentage_1000 = 100*(melt_1000 - melt_refer2)./melt_refer2;
display(length(percentage_1500));

percentage_1000 = [0, 0, 0, 0, 0, percentage_1000];
percentage_1500 = [0, 0, 0, 0, 0, percentage_1500];
%percentage_300 = [0, 0, 0, 0, 0, percentage_300];
%percentage_600 = [0, 0, 0, 0, 0, percentage_600];
time2 = linspace(0, 320, length(percentage_1500));
display(time2);

x = (T_1000.time_ky -T_1000.time_ky(1))/1e3;
xi = linspace(min(x), max(x), 20);  
yi = interp1(x, T_1000.percentage, xi, 'spline', 'extrap');

name = char(figname(end));

%% plot figure 
figure();

pos = get(gcf,'Position');
width = pos(3) * 3;
height = pos(4) * 1.5;
set(gcf,'Position',[pos(1) pos(2) width height]);
yyaxis left
if volume == 1 % plot melting volume 
    bar(time2, melt_1500 - melt_refer, 'FaceColor', [0.8500 0.3250 0.0980]);
    hold on;
    bar(time2, melt_1000 - melt_refer, 'FaceColor', [0.8500 0.5250 0.0980]);
    bar(time2, melt_600 - melt_refer, 'FaceColor', [0.8500 0.7250 0.0980]);
    bar(time2, melt_300 - melt_refer, 'FaceColor', [0.8500 0.9250 0.0980]);
%    ylim([120, 200]);
%     plot ((T_1000.time_ky -T_1000.time_ky(1))/1e3, T_1000.melt_v - T_1000.melt_refer,'linewidth',2);
%     hold on;
%     plot ((T_600.time_ky -T_600.time_ky(1))/1e3, T_600.melt_v - T_600.melt_refer,'linewidth',2);
%     plot ((T_300.time_ky -T_300.time_ky(1))/1e3, T_300.melt_v - T_300.melt_refer,'linewidth',2);
    ylabel({'Difference in extracted',...
    'melt  (km^3/km)'});
    output_image = ['output/AfricaModels2022/Figure/extraction/',name(1:end-4),...
     '_melting volume',];
elseif volume == 2  % plot melting percentage
    bar(time2, -0.5* percentage_1500, 'FaceColor', [0.8500 0.3250 0.0980]);
    hold on;
    bar(time2, percentage_1000, 'FaceColor', [0.8500 0.5250 0.0980]);
    bar(time1, percentage_600, 'FaceColor', [0.8500 0.7250 0.0980]);
    bar(time1, percentage_300, 'FaceColor', [0.8500 0.9250 0.0980]);
%    plot(time2, -0.5* percentage_1500, 'Color', [0.8500 0.3250 0.0980]);
%    hold on;
%    plot(time2, percentage_1000, 'Color', [0.8500 0.5250 0.0980]);
%    plot(time1, percentage_600, 'Color', [0.8500 0.7250 0.6980]);
%    plot(time1, percentage_300, 'Color', [0.8500 0.9250 0.3980]);
%    plot(time2, -0.3*percentage_600, 'Color', [0.8500 0.7250 0.0980]);
%    plot(time2, -0.1*percentage_300, 'Color', [0.8500 0.9250 0.0980]);

    ylabel({'Difference in extracted melt (%)'});
    output_image = ['output/AfricaModels2022/Figure/extraction/',name(1:end-4),...
     '_melting percentage_factor',];
else
    disp('please choose plot melting volume or percentage')
end 
% difference
% plot (T.time_ky/1e3, T.melt_v - T.melt_refer,'linewidth',2);
% difference in percentage 
% plot (T.time_ky/1e3, T.percentage*100,'linewidth',2);

xlabel('Model time (kyr)');
%ylabel('Melt volume, (km^2/km)');
% xlim([0.05 0.27])
% ylim([-1 3])

%water_range = max(T_1000.water_mass_kg) - min(T_1000.water_mass_kg);
yyaxis right
plot (T_1000.time_ky - T_1000.time_ky(1), 100*T_1000.water_mass_kg/T_1000.water_mass_kg(1),'linewidth',4, 'color', [0 0.4470 0.7410]);
ylabel('Change in Water mass (%)');
% ylabel('Water mass change (%)');
% legend({'differential melt with lake level drop',...
%     'surface water mass'}, 'Location','southeast');
legend({'max lake level change-1200m','max lake level change-800m',...
    'max lake level change-400m', 'max lake level change-200m'}, 'Location','northwest');

set(gca,'FontSize',25);
ax = gca;
ax.YAxis(1).Color = [0.8500 0.3250 0.0980];
ax.YAxis(2).Color = [0 0.4470 0.7410];
hold off;


%% export 
% % Create a png

print('-dpng','-r200',[output_image,'_factor.png']);

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
