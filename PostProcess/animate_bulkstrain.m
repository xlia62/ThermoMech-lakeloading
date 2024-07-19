function [mx,my,markstrain]=animate_bulkstrain(inpath,nfreq,ntotal,npx,npy,cmap,tempc,xlims,ylims)

% Inpath = folder where grid_###.mat and markers_###.mat files exist
% nfreq  = frequency of time stpes when files are saved
% ntotal = total number of time steps
% npx,npy = number of image pixels in x and y direction
% cmap = colormap
% tempc = temperature contours, ommit adding contours if empty

if nfreq == 1
    k = [1 nfreq:nfreq:ntotal];
else
    k = [nfreq:nfreq:ntotal];
end

yr2sec = 365.25*24*3600;

for j = k
    
    % Load files
    infile = [inpath,'/markers_',num2str(j),'.mat'];
    load(infile,'MX','MY','MBII','MI');
    infile = [inpath,'/grids_',num2str(j),'.mat'];
    load(infile);
    
    % Replace Air or Water with NaN
    i = MI == 1 | MI == 2;
    MBII(i) = NaN;
    
    if ~isempty(xlims)
       i = MX/1000 >= xlims(1) & MX/1000 <= xlims(2);
       ox = xlims(1);
       mx = MX(i) - ox*1000;
       my = MY(i);
       mbii = MBII(i);
    else
       xlims = [G.gridx(1) G.gridx(end)]/1000;
       ox = 0;
       mx = MX;
       my = MY;
       mbii = MBII;
    end
    if ~isempty(ylims)
        i = my/1000 >= ylims(1) & my/1000 <= ylims(2);
        oy = ylims(1);
        my = my(i) - oy*1000;
        mx = mx(i);
        mbii = mbii(i);
    else
        ylims = [G.gridy(1) G.gridy(end)]/1000;
        oy = 0;
        my = my;
    end
    
    
    [mx,my,markstrain] = vis_markers((xlims-ox)*1000,(ylims-oy)*1000,npx,npy,mx,my,mbii);
    
    mx = mx + ox*1000;
    my = my + oy*1000;
    
    if j == 1
        xlimits = [mx(1) mx(end)];
        ylimits = [my(1) my(end)];
    end
    
    % Plot the material proporties
    clf
    pcolor(mx/1000,my/1000,-log10(markstrain));
    hold on
    if ~isempty(tempc)
        [C,h]=contour(G.gridx/1000,G.gridy/1000,G.tk1-273,tempc,'k','linewidth',1);
    end
    colormap(cmap);
    shading flat;
%    set(gca,'renderd
    
    set(gca,'Ydir','reverse')
    %set(gca,'Clim',[0 1])
    title(['Bulk strain and Temperature, Step=',num2str(j),' Myr=',num2str(timesum*1e-6/(yr2sec))]);
    xlabel('x, km')
    ylabel('y, km')
    %set(gca,'clim',[-3 0])
    drawnow
end

if ~isempty(tempc)
    clabel(C,h,'manual','rotation',0)
end

end

