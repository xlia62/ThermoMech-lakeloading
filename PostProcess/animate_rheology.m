function animate_rheology(inpath,nfreq,ntotal,npx,npy,cmap,tempc,xlims,ylims)

% Inpath = folder where grid_###.mat and markers_###.mat files exist
% nfreq  = frequency of time stpes when files are saved
% ntotal = total number of time steps
% npx,npy = number of image pixels in x and y direction
% cmap = colormap
% tempc = temperature contours, ommit adding contours if empty

if nfreq > 1
    k = [1 nfreq:nfreq:ntotal];
else
    k = [nfreq:nfreq:ntotal];
end

yr2sec = 365.25*24*3600;

figure
% pos = get(gcf,'Position');
% width = pos(3) * 1.5;
% height = pos(4) * 1.1;
% set(gcf,'Position',[pos(1) pos(2) width height]);

%3000 x 600
% pos = get(gcf,'Position');
% width = pos(3) * 3;
% height = pos(4) * 1.5;
% set(gcf,'Position',[pos(1) pos(2) width height]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Optional: write to video
% video_name = [inpath,'/Rheology_', inpath];
% video = VideoWriter(video_name);
% try
%     close(video);
% catch
%     'here'
% end
% video.FrameRate = 10;   % if every timestep is saved, a little fast
% open(video);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for j = k
    
    % Load files
    infile = [inpath,'/markers_',num2str(j),'.mat'];
    load(infile,'MX','MY','META','MI');
    infile = [inpath,'/grids_',num2str(j),'.mat'];
    load(infile);
    
    % Replace Air or Water with NaN
    i = MI == 1 | MI == 2;
    META(i) = NaN;
    
    if ~isempty(xlims)
       i = MX/1000 >= xlims(1) & MX/1000 <= xlims(2);
       ox = xlims(1);
       mx = MX(i) - ox*1000;
       my = MY(i);
       meta = META(i);
    else
       xlims = [G.gridx(1) G.gridx(end)]/1000;
       ox = 0;
       mx = MX;
       my = MY;
       meta = META;
    end
    if ~isempty(ylims)
        i = my/1000 >= ylims(1) & my/1000 <= ylims(2);
        oy = ylims(1);
        my = my(i) - oy*1000;
        mx = mx(i);
        meta = meta(i);
    else
        ylims = [G.gridy(1) G.gridy(end)]/1000;
        oy = 0;
        my = my;
    end
    
    
    [mx,my,markrheology] = vis_markers((xlims-ox)*1000,(ylims-oy)*1000,npx,npy,mx,my,meta);
    
    mx = mx + ox*1000;
    my = my + oy*1000;
    
    if j == 1
        xlimits = [mx(1) mx(end)];
        ylimits = [my(1) my(end)];
    end
    
    % Plot the material proporties
    clf
    pcolor(mx/1000,my/1000,log10(markrheology));
    hold on
    if ~isempty(tempc)
        [C,h]=contour(G.gridx/1000,G.gridy/1000,G.tk1-273,tempc,'k','linewidth',1);
    end
     colors1 = [.7 .7 .7; .4 .4 .4; .005 .005 .005];
    n = 100;
    m = size(colors1,1);
    t0 = linspace(0,1,m)';
    t = linspace(0,1,n)';
    r = interp1(t0,colors1(:,1),t);
    g = interp1(t0,colors1(:,2),t);
    b = interp1(t0,colors1(:,3),t);
    cmap2 = [r,g,b];
    colormap(cmap);
    shading interp;
    
    set(gca,'Ydir','reverse','FontSize',18,'color',[0.9 0.9 0.9])
    set(gcf,'color',[0.9 0.9 0.9])
    set(gca,'Clim',[17 25])
    title(['Rheology, Step=',num2str(j),' Myr=',num2str(timesum*1e-6/(yr2sec))],'FontSize', 24);
    xlabel('Horizontal Position (km)','FontSize', 20,'FontWeight','bold')
    ylabel('Depth (km)', 'FontSize',20,'FontWeight','bold')
    c = colorbar;
    ylabel(c,'Viscosity (10^{{\itx}} Pa s)','FontSize',20);
    axis image
        if ~isempty(xlims)
        xlim(xlims)
    end
    if ~isempty(ylims)
        ylim(ylims)
    end
    drawnow
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % For video
%     pause(.2)
%     f = getframe(gcf);
%     video.writeVideo(f);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end

if ~isempty(tempc)
    clabel(C,h,'manual','rotation',0)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for video
% close(video)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

