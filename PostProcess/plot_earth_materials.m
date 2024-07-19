function output_image=plot_t_meltV(inpath)
% this is function plotting the temperature and strength

% input:
% Inpath = folder where grid_###.mat and markers_###.mat files exist

% output: 
% image of earth materials 

% required:
% NAN

% main code:
%% load data and define indexes

% Load files
infile = [inpath,'/litho_strength.mat'];
load(infile);
 

%% plot figure 
%figure();
ax_position = [0.15, 0.175, 0.73, 0.625];
figure('Position', [10 10 800 1000]);
% Create the first axes
hax1 = axes('YDir','reverse', 'FontSize',25);

% Plot something here

hplot1 = line(-1e-9*SE.sigd_tension, SE.z/1e3, 'LineWidth',4);
    
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

hplot2 = line(SE.T- 273, SE.z/1e3, 'LineWidth',4, 'Color', 'r', 'Parent', hax2);

% Link the y limits and position together
linkprop([hax1, hax2], {'ylim', 'Position'});
set(gca, 'Position', ax_position)

% Draw some labels
xlabel(hax1, 'Strengh (GPa)')
xlabel(hax2, ['Temperature (', char(176),'C)'])
ylabel(hax1, 'Depth (km)')

% Add a legend? Why not?!
legend([hplot1, hplot2], {'Strengh', 'Temperature'})

% set the size of figure

% set(gca, 'Position', [10 10 500 800])
% set(gca,'OuterPosition',[10 10 500 800]);
%% export 
figname = split(inpath, '/');
output_image = ['output/AfricaModels2022/Figure/extraction/',char(figname(end)),'_strength3'];
print('-dpng','-r200',[output_image,'.png'])


end
