function animate_matprop(inpath,nfreq,ntotal,npx,npy,cmap,tempc,xlims,ylims,movie_out)

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

%k([45 49 53 57 58]) = [];

yr2sec = 365.25*24*3600;

figure('Name',['RockProperties:', inpath],'NumberTitle','off')

% Optional: write to video
if movie_out == 1
    video_name = [inpath,'/Rocks'];
    video = VideoWriter(video_name,'MPEG-4');
    video.FrameRate = 5;   % if every timestep is saved, a little fast
    open(video);
end

% Backward compatibility, unwrap colormap if it is a cell
if iscell(cmap)
   cb_labels = cmap{2};
   cmap = cmap{1};
else
   cb_labels = {'Air','Water','Sed.','UOC','LOC','UCC','LCC','OLM','CLM','AsthM','HydrM','UOP','LOP','Plume','Melt','UM'};
end

count = 0;
ncount = length(k);

for j = k
    count = count + 1;
    % Load files
    infile = [inpath,'/markers_',num2str(j),'.mat'];
    
    chk = whos('-file',infile);
    
    load(infile,'MX','MY','MI');
    
    % Backward compatibility: Melt fraction
    if ismember('MXM',{chk.name})
        load(infile,'MXM');
    else
        MXM = zeros(size(MI));
    end
    
    % Backward compatibility: Fraction of eclogite
    if ismember('MECL',{chk.name})
        load(infile,'MECL');
    else
        MECL = zeros(size(MI));
    end
    
    infile = [inpath,'/grids_',num2str(j),'.mat'];
    load(infile);
    
    % Crop to desired area
    if ~isempty(xlims)
        i = MX/1000 >= xlims(1) & MX/1000 <= xlims(2);
        ox = xlims(1);
        mx = MX(i) - ox*1000;
        my = MY(i);
        mi = MI(i);
        mxm = MXM(i);
        mecl = MECL(i,1);  % MECL is a marknum x 5 array
    else
        xlims = [G.gridx(1) G.gridx(end)]/1000;
        ox = 0;
        mx = MX;
        my = MY;
        mi = MI;
        mxm = MXM;
        mecl = MECL(:,1); % MECL is a marknum x 5 array
    end
    
    if ~isempty(ylims)
        i = my/1000 >= ylims(1) & my/1000 <= ylims(2);
        oy = ylims(1);
        my = my(i) - oy*1000;
        mx = mx(i);
        mi = mi(i);
        mxm = mxm(i);
        mecl = mecl(i);
    else
        ylims = [G.gridy(1) G.gridy(end)]/1000;
        oy = 0;
    end
    
    % Change marker type to 15 if melt exist
    k = mxm > 0;
    mi(k) = 15;
    
    % Change marker type to 17 if eclogite exist
    k = mecl > 0;
    mi(k) = 17;
    
    % Change marker type to 18 if eclogite > 50%
    k = mecl > 0.5;
    mi(k) = 18;
    
    [mx,my,markcom] = vis_markers((xlims-ox)*1000,(ylims-oy)*1000,npx,npy,mx,my,mi);
    
    mx = mx + ox*1000;
    my = my + oy*1000;
    
    if j == 1
        xlimits = [mx(1) mx(end)];
        ylimits = [my(1) my(end)];
    end
    
    % Plot the material proporties   
    
    clf
    pcolor(mx/1000,my/1000,markcom);
    hold on
    if ~isempty(tempc)
        [C,h]=contour(G.gridx/1000,G.gridy/1000,G.tk1-273,tempc,'k','linewidth',1);
    end
    colormap(cmap);
    shading flat;
       
         %clabel(C,h,'LabelSpacing',600,'FontSize',12,'Margin',6);
%    clabel(C,h,'manual','FontSize',12,'Margin',6)
    
    % Plot title and labels
    set(gca,'Ydir','reverse','FontSize',12,'color',[0.9 0.9 0.9])
    set(gcf,'color',[0.9 0.9 0.9])
    set(gca,'Clim',[1 length(cmap)])
    
    timesum = timesum/yr2sec;
    if timesum>=1e6
        txt_time = [num2str(round(timesum*1e-6,2)),' Myr'];
    elseif timesum>=1e3
        txt_time = [num2str(round(timesum*1e-3,2)),' kyr'];
    else
        txt_time = [num2str(round(timesum,2)),' yr'];
    end
    
    title(['Material Properties: ',txt_time],'FontSize', 16);
    xlabel('Horizontal Position (km)','FontSize', 12,'FontWeight','bold')
    ylabel('Depth (km)', 'FontSize',12,'FontWeight','bold') 
    
    % Assemble the colorbar and its labels
    incr = (max(length(cmap))-1)/max(length(cmap));
%     colorbar('Ticks',1+incr/2:incr:max(length(cmap)),...
%         'TickLabels',cb_labels,...
%         'FontSize',12,'ydir','reverse');
    
   % axis image
    
    if ~isempty(xlims)
        xlim(xlims)
    end
    if ~isempty(ylims)
        ylim(ylims)
    end
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

