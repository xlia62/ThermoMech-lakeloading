function output_image=plot_PT_curve(inpath, step, zoomin)
% this is function plotting the temperature and strength

% input:
% Inpath = folder where grid_###.mat and markers_###.mat files exist
% zoomin = 0 or 1

% output: 
% image of earth materials 

% required:
% NAN

% main code:
%% load data and define indexes

% Load files
% infile = [inpath,'/litho_strength.mat'];
% load(infile);
infile = [inpath,'/grids_',num2str(step),'.mat']; % as G
load(infile);
infile = [inpath,'/markers_',num2str(step),'.mat'];
load(infile);
% load reference data when the lake doesn't change
figname = split(inpath, '/'); 
refername = replace(char(figname(end)), 'drop', 'nochange'); 
ref= load (['output/AfricaModels2022/Lake_model_extraction/', refername, '/grids_',num2str(step),'.mat']);
% zoomin = 1;
yr2sec = 365.25*24*3600;

% grid used in the function;
% gridx: [271×1 double]
% gridy: [221×1 double]
% MGII1 = G.rho1;
% unit = 'Density (km/m^3)';
% MGII1 = G.pr1;
% unit = 'Pressure (Pa)';
% MGII1 = G.tk1;
% unit = 'Temperature (K)';

%% plot figure 

% find the x index of max melting 
melt_max = max(MXM);
% id of max melt 
melt_max_i = find (MXM==melt_max);
% x of max melt
melt_x = MX(melt_max_i); % m
% find the grid of melt_x
tolerance = 500; %m 
xi = find(G.gridx > melt_x - tolerance & G.gridx < melt_x + tolerance);
% just use the first one for now
xi = xi(1);
display(melt_x/1e3, ' km');

figure('Position', [10 10 700 900]);
% plot P-T diagram
temperature = G.tk1(1:end-1,xi); % in K
pressure = G.pr1(:, xi); % in Pa
line (temperature -273, pressure/1e9, 'Color', [0 0.4470 0.7410],'LineWidth',2);
hold on; 
line (ref.G.tk1(1:end-1,xi) -273, ref.G.pr1(:, xi)/1e9,'LineWidth',2, 'Color', [0.4660 0.6740 0.1880]);

% plot melting curve 
%melting function from Melt_fraction
% 8,9 = Lithospheric mantle (dry): latent heat 400 kJ/kg
% case {8,9}
%     % Solidus Temperature
%     if (P<10000)
%         ts=1394+0.132899*P-0.000005104*P^2;
%     else
%         ts=2212+0.030819*(P-10000);
%     end
%     % Liquidus temperature
%     tl=2073+0.114*P;
%     % Latent heat
%     HL=400000;
% 
% % 11 = Hydrated mantle (wet): latent heat 400 kJ/kg
% case 11
%     % Solidus Temperature
%     if (P<2400)
%         ts=1240+49800/(P+323);
%     else
%         ts=1266-0.0118*P+0.0000035*P^2;
%     end
%     % Liquidus temperature
%     tl=2073+0.114*P;
%     % Latent heat
%     HL=400000;
            
% melt pressure 
melt_p = linspace(min(pressure), max(pressure), 200); % in pa
P = melt_p*1e-6; % MPa
% for dry olivine mantle 
solidus_t = 1394+0.132899*P-0.000005104*P.^2 -273; % in degree C 
liquidus_t = 1780 + 0.045*P - 2e-6*P.^2; % in degree C

line (solidus_t, melt_p/1e9, 'Color', 'red','LineWidth',2); 
line (liquidus_t, melt_p/1e9, 'Color', 'red', 'LineStyle','--','LineWidth',2); 

set(gca,'Ydir','reverse', 'FontSize',25)
% title(['Step:',num2str(step),' Model time:',num2str(timesum*1e-6/(yr2sec)), ' Myr']);
xlabel(['Temperature (', char(176),'C)'])
ylabel('Pressure (GPa)')
xlim([0. 2200]);
ylim([0. 4]);
if zoomin ~= 1
    legend( {'lake level drop', 'no change of lake level'},'Location','southwest')
end
if zoomin == 1
    xlim([1200 1400]);
    ylim([0.8 1.3]);
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure('Position', [10 10 700 900]);
% % Create the first axes
% hax1 = axes('YDir','reverse', 'FontSize',20);
% 
% % Plot something here
% 
% hplot1 = line(SE.sigd, SE.z/1e3);
%     
% % Create a transparent axes on top of the first one with it's xaxis on top
% % and no ytick marks (or labels)
% hax2 = axes('Position', get(hax1, 'Position'), ...  % Copy position
%             'XAxisLocation', 'top', ...             % Put the x axis on top
%             'YAxisLocation', 'right', ...           % Doesn't really matter
%             'Color', 'none', ...                    % Make it transparent
%             'YTick', [], ...                        % Don't show markers on y axis
%             'YDir','reverse',...
%             'FontSize',18);                          
% 
%         
% % Plot data with a different x-range here
% 
% hplot2 = line(SE.T- 273, SE.z/1e3, 'Color', 'r', 'Parent', hax2);
% 
% % Link the y limits and position together
% linkprop([hax1, hax2], {'ylim', 'Position'});
% 
% % Draw some labels
% xlabel(hax1, 'Strengh (Pa)')
% xlabel(hax2, ['Temperature (', char(176),'C)'])
% ylabel(hax1, 'Depth (km)')
% 
% % Add a legend? Why not?!
% legend([hplot1, hplot2], {'Strengh', 'Temperature'})
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% export 

output_image = ['output/AfricaModels2022/Figure/extraction/',char(figname(end)),'_PT_plot_',num2str(step)];
if zoomin == 1
    output_image = ['output/AfricaModels2022/Figure/extraction/',char(figname(end)),'_PT_plot_',num2str(step), '_zoomin(zoomin)'];
end
print('-dpng','-r200',[output_image,'.eps'])


end
