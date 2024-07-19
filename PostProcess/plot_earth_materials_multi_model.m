function output_image=plot_t_meltV(inpath1, inpath2)
% this is function plotting the temperature and strength

% input:
% Inpath1 = folder where grid_###.mat and markers_###.mat files exist
% Inpath2 = folder where grid_###.mat and markers_###.mat files exist
% type:
% 1: strength  2: temperature 

% output: 
% image of earth materials 

% required:
% NAN

% main code:
%% load data and define indexes

% Load files
infile = [inpath1,'/litho_strength.mat'];
a = load(infile);
infile = [inpath2,'/litho_strength.mat'];
b = load(infile); 

%% plot figure 
%figure();
ax_position = [0.15, 0.175, 0.73, 0.625];
figure('Position', [10 10 800 1000]);
% Create the first axes
hax1 = axes('YDir','reverse', 'FontSize',25);

% Plot something here

hplot1 = line(-1e-9*a.SE.sigd_tension, a.SE.z/1e3, 'LineWidth',4, 'Color',[0 0.4470 0.7410]);
hplot3 = line(-1e-9*b.SE.sigd_tension, b.SE.z/1e3, 'LineWidth',4, 'Color',[0 0.4070 0.3410]);    
% Create a transparent axes on top of the first one with it's xaxis on top
% and no ytick marks (or labels)
hax2 = axes('Position', get(hax1, 'Position'), ...  % Copy position
            'XAxisLocation', 'top', ...             % Put the x axis on top
            'YAxisLocation', 'right', ...           % Doesn't really matter
            'Color', 'none', ...                    % Make it transparent
            'YTick', [], ...                        % Don't show markers on y axis
            'YDir','reverse',...
            'FontSize',25);                          

set(gca, 'Position', ax_position)
% Plot data with a different x-range here

hplot2 = line(a.SE.T- 273, a.SE.z/1e3, 'LineWidth',4, 'Color', 'r', 'Parent', hax2);
hplot4 = line(b.SE.T- 273, b.SE.z/1e3, 'LineWidth',4, 'Color',  [0.8500 0.3250 0.0980], 'Parent', hax2);
% Link the y limits and position together
linkprop([hax1, hax2], {'ylim', 'Position'});
set(gca, 'Position', ax_position)

% Draw some labels
xlabel(hax1, 'Strengh (GPa)')
xlabel(hax2, ['Temperature (', char(176),'C)'])
ylabel(hax1, 'Depth (km)')

% Add a legend? Why not?!
%legend([hplot1, hplot2], {'Strengh', 'Temperature'})
%legend([hplot1, hplot3, hplot2, hplot4], {'Strengh (TL = 100km)', ...
%    'Strengh (TL  = 105km)', 'Temperature (TL = 100km)',...
%    'Temperature(TL= 105km)'},'Location','southwest', 'FontSize',22)
% 
% legend([hplot1, hplot3, hplot2, hplot4], {['Strengh (Tp = 1343'  char(176),'C)'], ...
%     ['Strengh (Tp = 1353'  char(176),'C)']', ['Temperature (Tp = 1343'  char(176),'C)']',...
%     ['Temperature (Tp = 1353'  char(176),'C)']},'Location','southwest', 'FontSize',22)

% legend([hplot1, hplot3, hplot2, hplot4], {'Strengh (reference crust)', ...
%     'Strengh (strong crust)', 'Temperature (reference crust)',...
%     'Temperature(strong crust)'},'Location','southwest', 'FontSize',22)

 legend([hplot1, hplot3, hplot2, hplot4], {'Strengh (reference crust)', ...
     'Strengh (weak crust)', 'Temperature (reference crust)',...
     'Temperature(weak crust)'},'Location','southwest', 'FontSize',22)

set(gca,'FontSize',21);
% set the size of figure

% set(gca, 'Position', [10 10 500 800])
% set(gca,'OuterPosition',[10 10 500 800]);
%% export 
figname1 = split(inpath1, '/');
figname2 = split(inpath2, '/');
output_image = ['output/AfricaModels2022/Figure/extraction/strength_comparison_',char(figname1(end)),...
    '_vs_', char(figname2(end))];
print('-dpng','-r200',[output_image,'.png'])


end
