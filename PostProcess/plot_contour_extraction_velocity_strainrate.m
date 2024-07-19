function animate_plasticstrain(inpath,step,xlims,ylims, MPROP )
%% 

% Inpath = folder where grid_###.mat and markers_###.mat files exist

% nfreq  = frequency of time stpes when files are saved
% ntotal = total number of time steps
% npx,npy = number of image pixels in x and y direction
% cmap = colormap
% tempc = temperature contours, ommit adding contours if empty
% inpath = 'output/AfricaModels2022/fastscape_line2';


npx = 10*(xlims(2) - xlims(1)); % xlims in km, so pixel will have 100 m by 100m resolution  
npy = 10*(ylims(2) - ylims(1));


j = step;
cmap = jet; %colormap(parula(5))

% if nfreq == 1
%    k = [1 nfreq:nfreq:ntotal];
% else
%     k = [nfreq:nfreq:ntotal];
% end



yr2sec = 365.25*24*3600;

    
% Load files
infile = [inpath,'/markers_',num2str(j),'.mat'];
load(infile);
infile = [inpath,'/grids_',num2str(j),'.mat']; % as G
load(infile);


 %% strain in markers
switch MPROP
    case 'META' % viscosity, Pa s
        MGII = META;
        unit = 'viscosity (Pa s)';
    case 'MGII' % Accumulated plastic strain
        MGII = MGII;
        unit = 'Accumulated plastic strain';
    case 'MBII' % Accumulated bulk strain
        MGII = MBII;
        unit = 'Accumulated bulk strain';
    case 'MPR' % pressure
        MGII = MPR;
        unit = 'Pressure (Pa)';
    case 'MTK' % Temperature, K
        MGII = MTK;
        unit = 'Temperature (K)';
    case 'MXM' % Cummulative Melt Fraction   
        MGII = MXM;
        unit = 'Cummulative Melt Fraction (%)';
    case 'MEXTC' % Cummulative Extracted Melt Fraction
        MGII = MEXTC;
        unit = {['Cummulative Extracted',newline,...
            'Melt Fraction (%)']};
    case 'MEXT' %  Extracted Melt Fraction
        MGII = MEXT;
        unit = 'Extracted Melt Fraction';
    case 'mextout' %  Extracted Melt Fraction
        MGII = mextout;
        unit = 'Extracted Melt Fraction';
    case 'MSXY' % Cummulative Extracted Melt Fraction
        MGII = MSXY;
        unit = 'shear stress (Pa)';
    case 'MSXX' % Cummulative Extracted Melt Fraction
        MGII = MSXX;
        unit = 'normal stress (Pa)';
    case 'MEXX' % normal strain rate
        MGII = MSXX;
        unit = 'normal strain rate (1/s)';
    case 'MEXY' % shear strain rate
        MGII = MSXX;
        unit = 'shear strain rate (1/s)';
    case 'MEII' % shear strain rate
        MGII = sqrt(MEXX.*MEXX + MEXY.*MEXY);
        unit = 'strain rate (1/s)';  
       
end


% Replace Air or Water or mantle with NaN
%i = MI == 1 | MI == 2| MI == 9| MI == 10;
i = MI == 1 | MI == 2;
MGII(i) = NaN;

% select based on x and y limits
if ~isempty(xlims)
    i = MX/1000 >= xlims(1) & MX/1000 <= xlims(2);
    ox = xlims(1);
    mx = MX(i) - ox*1000;
    my = MY(i);
    mgii = MGII(i);

else
    xlims = [G.gridx(1) G.gridx(end)]/1000;
    ox = 0;
    mx = MX;
    my = MY;
    mgii = MGII;

end
if ~isempty(ylims)
    i = my/1000 >= ylims(1) & my/1000 <= ylims(2);
    oy = ylims(1);
    my = my(i) - oy*1000;
    mx = mx(i);
    mgii = mgii(i);

else
    ylims = [G.gridy(1) G.gridy(end)]/1000;
    oy = 0;
    my = my;
end
    
% plot strain in grid 
[mx,my,markstrain] = vis_markers((xlims-ox)*1000,(ylims-oy)*1000,npx,npy,mx,my,mgii);

mx = mx + ox*1000;
my = my + oy*1000;
    
%% velocity in grids
vxb=zeros(G.ynum,G.xnum);
vyb=zeros(G.ynum,G.xnum);
% note that vx1 and vy1 are in m/s
[Xq,Yq] = meshgrid(linspace (xlims(1)*1e3, xlims(2)*1e3, npx),...
    linspace(ylims(1)*1e3, ylims(2)*1e3,npy));

for j=1:G.xnum
    for i=1:G.ynum
        vxb(i,j)=(G.vx1(i,j)+G.vx1(i+1,j))/2;
        vyb(i,j)=(G.vy1(i,j)+G.vy1(i,j+1))/2;
    end
end


Vx = interp2(G.gridx,G.gridy,vxb,Xq,Yq);
Vy = interp2(G.gridx,G.gridy,vyb,Xq,Yq);
Xq = Xq/1e3; % into km
Yq = Yq/1e3; % into km

%% Plot the material proporties
%clf
%figure('Position', [10 10 1000 800]);
% Position and size of figure on computer screen
plot_fig = 1;
if plot_fig
    
    figure();
    pos = get(gcf,'Position');
    width = pos(3) * 3;
    height = pos(4) * 1.5;
    set(gcf,'Position',[pos(1) pos(2) width height]);

    %pcolor(mx/1000,my/1000,-log10(abs(markstrain)))
    plot(linspace(min(mx)/1e3, max(mx)/1e3, 100), 1e-3*water_lev*ones(100), ...
        'color', 'blue', 'linewidth',4);
    %legend({'lake level'}, 'Location','northwest');
    hold on;
    % plot strain rate
    s = pcolor(mx/1000,my/1000, markstrain);
    colormap(cmap);
    c = colorbar; 
    s.FaceColor = 'interp';


    % [con, h] = contour(G.gridx/1000,G.gridy/1000,G.tk1-273,[500 800 1000 1200 1400 1600],'white','linewidth',1);
    % clabel(con, h,  'Fontsize', 22, 'Color', 'white', 'labelspacing', 400);
    %quiver(Xq(1:100:end),Yq(1:100:end),Vx(1:100:end),Vy(1:100:end),1,'k');



    %title(['Step=',num2str(j),' Model time: ',num2str(timesum*1e-6/(yr2sec)), 'Myr']);
    xlabel('Horizon (km)','FontSize',25)
    ylabel('Depth (km)','FontSize',25)
    ylabel(c,unit,'FontSize',25);
    shading interp;

    ax=gca;
    set(gca,'Ydir','reverse');
    % set up the scale for MGII

    ax.FontSize = 25;
    % if MPROP == 'MGII'
    %     caxis([0 10]);
    %end
    % if MPROP == 'MXM'
    %     caxis([0.14 0.16]);
    % end
    %caxis([8.7 9.7]);
    %caxis([1 12]);
    caxis([0, 1e-15]);
    %legend( {'lake level'}, 'FontSize',16)
    axis image
    if ~isempty(xlims)
        xlim(xlims)
    end
    if ~isempty(ylims)
        ylim(ylims)
    end
    drawnow


    % move the sticky air from figure
    move = gca; 
    tick = move.YTick - 20;
    yticklabels({tick});

    % Create a png
    figname = split(inpath, '/');
    output_image = ['output/AfricaModels2022/Figure/extraction/',char(figname(end)),'_step',num2str(step),'_velocity_', MPROP, '_zoom'];
    print('-dpng','-r200',[output_image,'.png'])
    
end

% export as mat 
figname = split(inpath, '/');
output_mat = ['output/AfricaModels2022/Mat/',char(figname(end)),'_step',num2str(step),...
    '_x_', num2str(xlims(1)), '-',  num2str(xlims(2)),...
    '_y_', num2str(ylims(1)), '-',  num2str(ylims(2)),'_', MPROP, '_zoom'];
% save ([output_mat,'.mat'], markstrain)
dlmwrite([output_mat,'.txt'], markstrain, '\t');
output_mat2 = ['output/AfricaModels2022/Mat/',char(figname(end)),'_step',num2str(step),...
    '_x_', num2str(xlims(1)), '-',  num2str(xlims(2)),...
    '_y_', num2str(ylims(1)), '-',  num2str(ylims(2)),'_velocity_y'];
dlmwrite([output_mat2,'.txt'], Vy, '\t');
output_mat3 = ['output/AfricaModels2022/Mat/',char(figname(end)),'_step',num2str(step),...
    '_x_', num2str(xlims(1)), '-',  num2str(xlims(2)),...
    '_y_', num2str(ylims(1)), '-',  num2str(ylims(2)),'_velocity_x'];
dlmwrite([output_mat3,'.txt'], Vx, '\t');


end

