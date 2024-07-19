function animate_temp(inpath,nfreq,ntotal,npx,npy,cmap,tempc,xlims,ylims,opt)

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
% Scaled for the 3500km x 1500km grid
figure
% pos = get(gcf,'Position');
% width = pos(3) * 1.5;
% height = pos(4) * 1.1;
% set(gcf,'Position',[pos(1) pos(2) width height]);

%3000 x 600
pos = get(gcf,'Position');
width = pos(3) * 3;
height = pos(4) * 1.5;
set(gcf,'Position',[pos(1) pos(2) width height]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Optional: write to video
% video_name = [inpath,'/Temperature_', inpath,'.mp4'];
% video = VideoWriter(video_name,'MPEG-4');
% try
%     close(video);
% catch
%     'here'
% end
% video.FrameRate = 10;   % if every timestep is saved, a little fast
% open(video);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[X,Y] = meshgrid(10:20:1190,20:10:550);

for j = k
    
    % Load files
    infile = [inpath,'/markers_',num2str(j),'.mat'];
    load(infile,'MX','MY','MTK','MI');
    infile = [inpath,'/grids_',num2str(j),'.mat'];
    load(infile);
   
    [Vx,Vy]=velocityplot(X*1000,Y*1000,G);
    MTK = MTK - 273;
    % Replace Air or Water with NaN
    i = MI == 1 | MI == 2;
    MTK(i) = NaN;
    
    if ~isempty(xlims)
       i = MX/1000 >= xlims(1) & MX/1000 <= xlims(2);
       ox = xlims(1);
       mx = MX(i) - ox*1000;
       my = MY(i);
       mtemp = MTK(i);
    else
       xlims = [G.gridx(1) G.gridx(end)]/1000;
       ox = 0;
       mx = MX;
       my = MY;
       mtemp = MTK;
    end
    if ~isempty(ylims)
        i = my/1000 >= ylims(1) & my/1000 <= ylims(2);
        oy = ylims(1);
        my = my(i) - oy*1000;
        mx = mx(i);
        mtemp = mtemp(i);
    else
        ylims = [G.gridy(1) G.gridy(end)]/1000;
        oy = 0;
        my = my;
    end
    
    
    [mx,my,markertemp] = vis_markers((xlims-ox)*1000,(ylims-oy)*1000,npx,npy,mx,my,mtemp);
    
    mx = mx + ox*1000;
    my = my + oy*1000;
    
    if j == 1
        xlimits = [mx(1) mx(end)];
        ylimits = [my(1) my(end)];
    end
    
    % Plot the material proporties
    clf
    if opt
        pcolor(mx/1000,my/1000,markertemp-repmat(markertemp(:,2),1,npx));
    else
        pcolor(mx/1000,my/1000,markertemp);
    end
    caxis([1000 1600])
    colormap(cmap);
    hold on
 %   curvvec(X,Y,Vx,Vy,'color','k','thin',3)
    
    % Optional - red-yellow intensity colormap
%     red_yellow = [0 0 0; 1 0 0; 1 193/255 37/255; 1 1 1];
%     n = 40;
%     m = size(red_yellow,1);
%     t0 = linspace(0,1,m)';
%     t = linspace(0,1,n)';
%     r = interp1(t0,red_yellow(:,1),t);
%     g = interp1(t0,red_yellow(:,2),t);
%     b = interp1(t0,red_yellow(:,3),t);
%     red_yellow = [r,g,b];
%     blue_cold = [0 0 1];
%     cmap2 = [blue_cold; red_yellow];
%     colormap(cmap2);
%     caxis([1800,2800])
    
    shading interp;

    
    set(gca,'Ydir','reverse','FontSize',18,'color',[0.9 0.9 0.9])
    set(gcf,'color',[0.9 0.9 0.9])
    title(['Temperature, Step=',num2str(j),' Myr=',num2str(timesum*1e-6/(yr2sec))],'FontSize', 24);
    xlabel('Horizontal Position (km)','FontSize', 20,'FontWeight','bold')
    ylabel('Depth (km)', 'FontSize',20,'FontWeight','bold')
    c = colorbar;
    ylabel(c,'Temperature (^{o}C)','FontSize',20);
    axis image
    drawnow
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % For video
%     pause(.2)
%     f = getframe(gcf);
%     video.writeVideo(f);
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for video
% close(video)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

