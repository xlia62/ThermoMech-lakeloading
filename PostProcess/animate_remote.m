function animate_remote(fig,cmap,tempc,movie_out,sticky_layer)

% Inpath = folder where grid_###.mat and markers_###.mat files exist
% nfreq  = frequency of time stpes when files are saved
% ntotal = total number of time steps
% npx,npy = number of image pixels in x and y direction
% cmap = colormap
% tempc = temperature contours, ommit adding contours if empty


figure('Name','RockProperties')
set(gcf,'Position',[250 75 1280 720])

% Optional: write to video
if movie_out == 1
    video_name = 'RocksProperties';
    video = VideoWriter(video_name,'MPEG-4');
    video.FrameRate = 4;   % if every timestep is saved, a little fast
    video.Quality = 100;
    open(video);
end

% Backward compatibility, unwrap colormap if it is a cell
if iscell(cmap)
   cb_labels = cmap{2};
   cmap = cmap{1};
else
   cb_labels = {'Air','Water','Sed.','UOC','LOC','UCC','LCC','OLM','CLM','AsthM','HydrM','UOP','LOP','Plume','Melt','UM'};
end

ncount = length(fig);

xlimits = [fig(1).mx(1) fig(1).mx(end)];
ylimits = [0 180];%[fig(1).my(1) fig(1).my(end)] - sticky_layer;

for j = 1:ncount    
    
    fig(j).my = fig(j).my - sticky_layer;
    
        % Plot the material proporties   
    
    clf
    pcolor(fig(j).mx,fig(j).my,fig(j).markcom);
    hold on
    if tempc && isfield(fig(j),'tempc')
            fig(j).gridy = fig(j).gridy - sticky_layer;
            [C,h]=contour(fig(j).gridx,fig(j).gridy,fig(j).tk1,fig(j).tempc,'color',[0 0.45 0.74],'linewidth',2);
   end
    colormap(cmap);
    shading flat;
       
%    clabel(C,h,'LabelSpacing',1500,'FontSize',20,'color',[0 0.45 0.74]);
    
    % Plot title and labels
    set(gca,'Ydir','reverse','FontSize',12,'color',[0.9 0.9 0.9])
    set(gcf,'color',[0.9 0.9 0.9])
    set(gca,'Clim',[1 length(cmap)])
    title(['Material Properties: ',num2str(round(fig(j).time,2)),' Myr'],'FontSize', 16);
    xlabel('Horizontal Position (km)','FontSize', 12,'FontWeight','bold')
    ylabel('Depth (km)', 'FontSize',12,'FontWeight','bold') 
    
        % tl =clabel(C,h,'manual','FontSize',20,'Margin',20);

    % Assemble the colorbar and its labels
    incr = (max(length(cmap))-1)/max(length(cmap));
%     colorbar('Ticks',1+incr/2:incr:max(length(cmap)),...
%         'TickLabels',cb_labels,...
%         'FontSize',12,'ydir','reverse');
    
   % axis image
    
    
  
  axis image
axis tight 
xlim(xlimits)   
   ylim(ylimits)
    drawnow
%             set(gcf,'Position',[118 253 1102 552])

    % For video
    if movie_out == 1
        
        pause(.2)
        f = getframe(gcf);
        video.writeVideo(f);
        
    elseif movie_out == 2
        set(gcf,'Position',[118 253 1102 552])
        f = getframe(gcf);
        if count == 1
            [im,map] = rgb2ind(f.cdata,256,'nodither');
            im(1,1,1,ncount) = 0;
        end
        im(:,:,1,count) = rgb2ind(f.cdata,map,'nodither');
                
    end
    
end

% if ~isempty(tempc) && ~movie_out
%     %clabel(C,h,'manual','rotation',0)
%     th = clabel(C,h,'manual','FontSize',12,'Margin',6);
%     set(th,'VerticalAlignment','bottom')
% elseif ~isempty(tempc)
%     clabel(C,h,'LabelSpacing',244,'FontSize',12,'Margin',6);
% end

% For video
if movie_out == 1
    close(video)
elseif movie_out == 2
    imwrite(im,map,'movie2.gif','DelayTime',0.2,'LoopCount',inf)
end


end

