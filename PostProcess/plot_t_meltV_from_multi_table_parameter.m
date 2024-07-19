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

factor = 0.8; %4.5 for lab65 for tp10 last half
%factor = 10; %3; % wet

% strong =  readtable('output/AfricaModels2022/Table/Lake_drop_restart_71_strongcrust_dry_melting_extract_history_112-174_dt_1_nomark_36_x_0-300_y_0-200.csv');
% viscous =  readtable('output/AfricaModels2022/Table/Lake_drop_restart_71_weakcrust_dry_melting_extract_history_504-572_dt_1_nomark_36_x_0-300_y_0-200.csv');
% tp10 =  readtable('output/AfricaModels2022/Table/Lake_drop_restart_71_1000m_Tp10_dry_melting_extract_history_530-598_dt_1_nomark_36_x_0-300_y_0-200_modify.csv');
% lab65 =  readtable('output/AfricaModels2022/Table/Lake_drop_restart_71_LAB65_dry_melting_extract_history_139-206_dt_1_nomark_36_x_0-300_y_0-200.csv');
% exten27 =  readtable('output/AfricaModels2022/Table/Lake_drop_restart_71_v27_dry_melting_extract_history_527-595_dt_1_nomark_36_x_0-300_y_0-200_modify.csv');
%highres = readtable('output/AfricaModels2022/Table/Lake_drop_restart_71_1000m_highres_melting_extract_history_620-688_dt_1_nomark_36_x_0-300_y_0-250.csv');
sed = readtable('output/AfricaModels2022/Table/Lake_drop_restart_71_dry_sed_melting_extract_history_197-259_dt_1_nomark_36_x_0-300_y_0-250.csv');

ref = readtable('output/AfricaModels2022/Table/Lake_drop_restart_71_dry_melting_extract_history_197-253_dt_1_nomark_36_x_0-300_y_0-200.csv');
figname = split('output/AfricaModels2022/Table/Lake_drop_restart_71_dry_melting_extract_history_197-253_dt_1_nomark_36_x_0-300_y_0-200', '/'); 


% consider if ploting every 2 element 
% note that the total number could be odd or even 
% melt_strong = mean(reshape(strong.melt_v(2:end), 2, []));
% melt_viscous = mean(reshape(viscous.melt_v(2:end), 2, []));
% melt_tp10 = mean(reshape(tp10.melt_v(2:end), 2, []));
% melt_lab65 = mean(reshape(lab65.melt_v(1:end), 2, []));
% melt_exten27 = mean(reshape(exten27.melt_v(2:end), 2, []));%
%melt_highres = mean(reshape(highres.melt_v(2:end), 2, []));
melt_sed = mean(reshape(sed.melt_v(2:end), 2, []));

% water mass value 
% water_strong = mean(reshape(strong.water_mass_kg(2:end), 2, []));
% water_viscous = mean(reshape(viscous.water_mass_kg(2:end), 2, []));
% water_tp10 = mean(reshape(tp10.water_mass_kg(2:end), 2, []));
% water_lab65 = mean(reshape(lab65.water_mass_kg(1:end), 2, []));
% water_exten27 = mean(reshape(exten27.water_mass_kg(2:end), 2, []));
%water_highres = mean(reshape(highres.water_mass_kg(2:end), 2, []));
water_sed = mean(reshape(sed.water_mass_kg(2:end), 2, []));

% the reference value for precentage calculation 
% melt_strong_r = mean(reshape(strong.melt_refer(2:end), 2, []));
% melt_viscous_r = mean(reshape(viscous.melt_refer(2:end), 2, []));
% melt_tp10_r = mean(reshape(tp10.melt_refer(2:end), 2, []));
% melt_lab65_r = mean(reshape(lab65.melt_refer(1:end), 2, []));
% melt_exten27_r = mean(reshape(exten27.melt_refer(2:end), 2, []));
%melt_highres_r = mean(reshape(highres.melt_refer(2:end), 2, []));
melt_sed_r = mean(reshape(sed.melt_refer(2:end), 2, []));


% the primary model 
melt_ref =  mean(reshape(ref.melt_v(2:end), 2, []));
melt_ref_r =  mean(reshape(ref.melt_refer(2:end), 2, []));

% time
% time_strong = mean(reshape(strong.time_ky(2:end), 2, []));
% time_viscous = mean(reshape(viscous.time_ky(2:end), 2, []));
% time_tp10 = mean(reshape(tp10.time_ky(2:end), 2, []));
% time_lab65 = mean(reshape(lab65.time_ky(1:end), 2, []));
% time_exten27 = mean(reshape(exten27.time_ky(2:end), 2, []));
time_ref = mean(reshape(ref.time_ky(2:end), 2, []));
%time_highres = mean(reshape(highres.time_ky(2:end), 2, []));
time_sed = mean(reshape(sed.time_ky(2:end), 2, []));

name = char(figname(end));

%% plot figure 
figure();

pos = get(gcf,'Position');
width = pos(3) * 3;
height = pos(4) * 1.5;
set(gcf,'Position',[pos(1) pos(2) width height]);
yyaxis left
if volume == 1 % plot melting volume 
    
    % tp10
%     bar(time_tp10 - time_tp10(1), (melt_tp10 - melt_tp10_r), 'FaceColor', [0.8500 0.3250 0.0980]);
%     hold on;
%     bar(time_ref - time_ref(1)-1, melt_ref - melt_ref_r, 'FaceColor', [0.8500 0.5250 0.0980]);
%     output_image = ['output/AfricaModels2022/Figure/extraction/',name(1:end-4),...
%      '_melting_volume_tp10_dry_modify',];
%  
%     %strong crust
%     bar(time_ref - time_ref(1), (melt_ref - melt_ref_r), 'FaceColor', [0.8500 0.5250 0.0980]);   
%     hold on;
%     bar(time_strong - time_strong(1), (melt_strong - melt_strong_r), 'FaceColor', [0.8500 0.3250 0.0980]);
%     output_image = ['output/AfricaModels2022/Figure/extraction/',name(1:end-4),...
%      '_melting_volume_strongcrust_dry',];
 
%     %viscous crust
%      bar(time_ref - time_ref(1), (melt_ref - melt_ref_r), 'FaceColor', [0.8500 0.5250 0.0980]);
%      hold on;
%      bar(time_viscous - time_viscous(1), (melt_viscous - melt_viscous_r), 'FaceColor', [0.8500 0.3250 0.0980]);  
%      output_image = ['output/AfricaModels2022/Figure/extraction/',name(1:end-4),...
%       '_melting_volume_weakcrust_dry',];

    %LAB depth   
%     bar(time_ref - time_ref(1), (melt_ref - melt_ref_r), 'FaceColor', [0.8500 0.5250 0.0980]);
%     hold on;
%     bar(time_lab65 - time_lab65(1), factor*(melt_lab65 - melt_lab65_r), 'FaceColor', [0.8500 0.3250 0.0980]);
%     output_image = ['output/AfricaModels2022/Figure/extraction/',name(1:end-4),...
%      '_melting_volume_lab65_dry',];

%     % extension rate  
%    bar(time_exten27 - time_exten27(1), melt_exten27 - melt_exten27_r, 'FaceColor', [0.8500 0.3250 0.0980]);
%    hold on;
%    bar(time_ref - time_ref(1), (melt_ref - melt_ref_r), 'FaceColor', [0.8500 0.5250 0.0980]);
%    output_image = ['output/AfricaModels2022/Figure/extraction/',name(1:end-4),...
%     '_melting_volume_exten27_dry_modify',];
 
     % extension rate after 71
%     bar(time_exten50 - time_exten50(1), melt_exten50 - melt_exten50_r, 'FaceColor', [0.8500 0.3250 0.0980]);
%     hold on;
%     bar(time_ref - time_ref(1), (melt_ref - melt_ref_r), 'FaceColor', [0.8500 0.5250 0.0980]);
%     output_image = ['output/AfricaModels2022/Figure/extraction/',name(1:end-4),...
%      '_melting_volume_exten50_after71',];
%  
    % high resolution 
%      bar(time_ref - time_ref(1), (melt_ref - melt_ref_r), 'FaceColor', [0.8500 0.5250 0.0980]);
%      hold on;
%      %%%%% NOOOOOOTTEEEE we have a factor here!!!!
%      bar(time_highres - time_highres(1), 0.5*(melt_highres- melt_highres_r), 'FaceColor', [0.8500 0.3250 0.0980]);  
%      output_image = ['output/AfricaModels2022/Figure/extraction/',name(1:end-4),...
%       '_melting_volume_highres_dry',];
% 
%     ylabel({'Difference in extracted',...
%     'melt (km^3/km)'});

    % sediment 
     bar(time_ref - time_ref(1), (melt_ref - melt_ref_r), 'FaceColor', [0.8500 0.5250 0.0980]);
     hold on;
     %%%%% 
     bar(time_sed - time_sed(1), (melt_sed- melt_sed_r), 'FaceColor', [0.8500 0.3250 0.0980]);  
     output_image = ['output/AfricaModels2022/Figure/extraction/',name(1:end-4),...
      '_melting_volume_sed_dry',];

    ylabel({'Difference in extracted',...
    'melt (km^3/km)'});

elseif volume == 2  % plot melting percentage
    % strong mantle
%     bar(time_strong - time_strong(1), factor*100*(melt_strong - melt_strong_r)./melt_strong_r, 'FaceColor', [0.8500 0.3250 0.0980]);
%     hold on;
%     bar(time_ref - time_ref(1), 100*(melt_ref - melt_ref_r)./melt_ref_r, 'FaceColor', [0.8500 0.5250 0.0980]);
    
    % extension
    %bar(time_exten0_1 - time_exten0_1(1), 100*(melt_exten0_1 - melt_exten0_1_r)./melt_exten0_1_r, 'FaceColor', [0.8500 0.7250 0.0980]);   
%     bar(time_exten27 - time_exten27(1), 100*(melt_exten27 - melt_exten27_r)./melt_exten27_r, 'FaceColor', [0.8500 0.3250 0.0980]);
%     hold on;
%     bar(time_ref - time_ref(1), 100*(melt_ref - melt_ref_r)./melt_ref_r, 'FaceColor', [0.8500 0.5250 0.0980]);
    
    %lab depth
%     bar(time_lab65 - time_lab65(1), 100*(melt_lab65 - melt_lab65_r)./melt_lab65_r, 'FaceColor', [0.8500 0.3250 0.0980]);
%     hold on;
%     bar(time_ref - time_ref(1), 100*(melt_ref - melt_ref_r)./melt_ref_r, 'FaceColor', [0.8500 0.5250 0.0980]);
%     

%     %tp100
%     bar(time_tp10 - time_tp10(1), factor*100*(melt_tp10 - melt_tp10_r)./melt_tp10_r, 'FaceColor', [0.8500 0.3250 0.0980]);
%     hold on;
%     bar(time_ref - time_ref(1), 100*(melt_ref - melt_ref_r)./melt_ref_r, 'FaceColor', [0.8500 0.5250 0.0980]);
% %     
%     
%     ylabel({'Differential extracted',...
%     'melting volume (%)'});
%     output_image = ['output/AfricaModels2022/Figure/extraction/',name(1:end-4),...
%      '_melting_percentage_tp10',];
%     % extension rate  
%     bar(time_exten27 - time_exten27(1), factor*100*(melt_exten27 - melt_exten27_r)./melt_exten27_r, 'FaceColor', [0.8500 0.3250 0.0980]);
%     hold on;
%     bar(time_ref - time_ref(1), factor*100*(melt_ref - melt_ref_r)./melt_ref_r, 'FaceColor', [0.8500 0.5250 0.0980]);
%     output_image = ['output/AfricaModels2022/Figure/extraction/',name(1:end-4),...
%      '_melting_percentage_exten27_dry',];


else
    disp('please choose plot melting volume or percentage')
end 
% difference
% plot (T.time_ky/1e3, T.melt_v - T.melt_refer,'linewidth',2);
% difference in percentage 
% plot (T.time_ky/1e3, T.percentage*100,'linewidth',2);
grid on;
xlabel('Model time (kyr)');
%ylabel('Melt volume, (km^3/km)');
% xlim([0.05 0.27])
%ylim([0 1.5])


yyaxis right
%plot ((ref.time_ky - ref.time_ky(1)), ref.water_mass_kg,'linewidth',4, 'color', [0 0.4470 0.7410]);
%plot (time_exten27 - time_exten27(1), water_exten27,'linewidth',4, 'color', [0 0.4470 0.7410]);
%plot (time_lab65 - time_lab65(1), water_lab65,'linewidth',4, 'color', [0 0.4470 0.7410]);
%plot (time_tp10 - time_tp10(1), water_tp10,'linewidth',4, 'color', [0 0.4470 0.7410]);
%plot (time_strong - time_strong(1), water_strong,'linewidth',4, 'color', [0 0.4470 0.7410]);
%plot (time_viscous - time_viscous(1), water_viscous,'linewidth',4, 'color', [0 0.4470 0.7410]);
%plot (time_sed - time_sed(1), water_sed,'linewidth',4, 'color', [0 0.4470 0.7410]);
plot (sed.time_ky - sed.time_ky(1), sed.water_mass_kg,'linewidth',4, 'color', [0 0.4470 0.7410]);

ylabel('Water mass (kg/km)');


% legend({'100% crust cohesion', '70% crust cohesion',...
%     'water level curve',  'water level curve2'}, 'Location','northwest');

% legend({'Tp = 1353 °C', 'Tp = 1343 °C',...
%     'water level curve', }, 'Location','northwest');
% legend({ 'reference crust','strong curst',...
%     'water level curve'}, 'Location','northwest');
%  legend({'reference crust','weak crust', ...
%      'water level curve'}, 'Location','northwest');
% legend({'inital lithosphere thickness = 100 km', 'inital lithosphere thickness = 105 km',...
%     'water level curve'}, 'Location','northwest');
%legend({'half extension rate= 3 mm/yr',...
%    'half extension rate= 2.5 mm/yr', 'water level curve'}, 'Location','northwest');
%  legend({'original resolution','fine resolution', ...
%      'water level curve'}, 'Location','northwest');
 legend({'deep lake','shallow lake', ...
      'water level curve'}, 'Location','northwest');

set(gca,'FontSize',25);
ax = gca;
ax.YAxis(1).Color = [0.8500 0.3250 0.0980];
ax.YAxis(2).Color = [0 0.4470 0.7410];

hold off;


%% export 
% % Create a png

print('-dpng','-r200',[output_image,'.png']);

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
