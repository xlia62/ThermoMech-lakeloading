function plot_t_meltV(inpath, melt_i, type)
% this is function plot the melting history(km^3/km)

% input:
% Inpath = folder where csv files exist
% melt_i choose to calculate melt or extract 1:MXM 2:MEXTC 

% type:  1: diff volume,  2: diff percentage, 3 plot two volume



% output: 
% image of differental melting history(km^3/km) abd 

% required:
% totalmelt_lake function

% main code:
%% load data and define indexes

T =  readtable(inpath);
% T_3000 =  readtable('output/AfricaModels2022/Table/Lake_drop_restart_90_2p_melting_history_300-369.csv');
figname = split(inpath, '/'); 
T.percentage = (T.melt_v - T.melt_refer)./(T.melt_refer) ; 

% consider if ploting every 2 element 
% note that the total number is odd
interval = 1;
start = 1; % just for 40 kyr period
melt_v2 = mean(reshape(T.melt_v(start:end), interval, []));
melt_refer2 = mean(reshape(T.melt_refer(start:end), interval, []));
time2 = mean(reshape((T.time_ky(start:end) - T.time_ky(start)), interval, []));
percentage2 = 100*(melt_v2 - melt_refer2)./melt_refer2;
%     
    
%% plot figure 
figure();

pos = get(gcf,'Position');
width = pos(3) * 3;
height = pos(4) * 1.5;
set(gcf,'Position',[pos(1) pos(2) width height]);

yyaxis left

% % time 
% time = (T.time_ky - T.time_ky(1))/1e3;
% % different melt V
% dif_v = T.melt_v - T.melt_refer;
% % smooth the V
% xx = linspace(min(time), max(time), 20);
% smooth_dif_v = spline (time, dif_v, xx);

if type == 1 % volume
%     plot ((T.time_ky - T.time_ky(1))/1e3, T.melt_v - T.melt_refer,'linewidth',2);
%     hold on;
  	%bar((T.time_ky - T.time_ky(1)), T.melt_v - T.melt_refer);
    bar(time2, melt_v2 - melt_refer2, 'FaceColor', [0.8500 0.3250 0.0980]);
elseif type == 2  % 
    %plot (T.time_ky - T.time_ky(1), T.percentage*100,'linewidth',2);
    %bar(time2, percentage2, 'FaceColor', [0.8500 0.3250 0.0980]);
    bar((T.time_ky - T.time_ky(1)), 100*T.percentage, 'FaceColor', [0.8500 0.3250 0.0980]); % this is just for 40 kyr period 
elseif type ==3
%     plot (T.time_ky - T.time_ky(1), T.melt_v,'linewidth',2);
%     hold on;
%     plot (T.time_ky - T.time_ky(1), T.melt_refer,'linewidth',2);
    bar(time2, melt_v2, 'FaceColor', [0.8500 0.3250 0.0980]);
    hold on;
    bar(time2, melt_refer2, 'FaceColor', [0.8500 0.5250 0.0980]);
    ylim([30, 60]);
else
    disp('wrong type')
end

xlabel('Model time (kyr)');
%ylabel('Melt volume, (km^2/km)');
if melt_i == 1 & type == 1
    ylabel('Difference in melt (km^3/km)');
elseif melt_i == 1 & type == 2
    ylabel('Difference in melt (%)');
elseif melt_i == 2 & type == 1
    ylabel({'Difference in extracted',...
    'melt (km^3/km)'});
elseif melt_i == 2 & type == 2
    ylabel({'Difference in ', ...
        'extracted melt (%)'});
elseif type == 3
    ylabel('Extracted melt (km^3/km)');
else 
    ylabel('wrong input');
end
% set(gca, 'YScale', 'log');
% xlim([0.05 0.27])
% ylim([8 25])

yyaxis right
plot ((T.time_ky - T.time_ky(1)), T.water_mass_kg,'linewidth',4, 'color', [0 0.4470 0.7410]);
ylabel('Water mass (kg/km)');

if type == 1 || type == 2 % volume
    legend({'difference in melt with lake level change',...
        'surface water mass'}, 'Location','southwest');
elseif type == 3
    legend({'melt with lake level change','melt with no lake level change',...
        'surface water volume'}, 'Location','southwest');
    
end
set(gca,'FontSize',25);
ax = gca;
ax.YAxis(1).Color = [0.8500 0.3250 0.0980];
ax.YAxis(2).Color = [0 0.4470 0.7410];


hold off;


%% export 
% % Create a png
name = char(figname(end));

output_image = ['output/AfricaModels2022/Figure/extraction/',name(1:end-4),...
     '_difference_', num2str(melt_i),'_type_test_',num2str(type)];
print('-dpng','-r200',[output_image,'.png']);



%print('-depsc','-r200',[output_image,'.eps']);
% 
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
