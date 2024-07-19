function animate_meltextc(inpath,nfreq,ntotal,npx,npy,cmap,tempc,xlims,ylims)

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

for j = k
    
    % Load files
    infile = [inpath,'/markers_',num2str(j),'.mat'];
    chk = whos('-file',infile);
    
    if ismember('MEXTC',{chk.name})
        load(infile,'MX','MY','MI','MTK','MEXTC');
    else
        load(infile,'MX','MY','MTK','MI');
        MEXTC = zeros(size(MI));
    end
    
    infile = [inpath,'/grids_',num2str(j),'.mat'];
    load(infile);
    
    % Replace Air or Water with NaN
    i = MI == 1 | MI == 2;
    MTK(i) = NaN;
    
    if ~isempty(xlims)
        i = MX/1000 >= xlims(1) & MX/1000 <= xlims(2);
       ox = xlims(1);
       mx = MX(i) - ox*1000;
       my = MY(i);
       mtemp = MTK(i);
       mextc = MEXTC(i);
    else
       xlims = [G.gridx(1) G.gridx(end)]/1000;
       ox = 0;
       mx = MX;
       my = MY;
       mtemp = MTK;
       mextc = MEXTC;

    end
    if ~isempty(ylims)
        i = my/1000 >= ylims(1) & my/1000 <= ylims(2);
        oy = ylims(1);
        my = my(i) - oy*1000;
        mx = mx(i);
        mtemp = mtemp(i);
        mextc = mextc(i);
    else
        ylims = [G.gridy(1) G.gridy(end)]/1000;
        oy = 0;
    end
    
    
    [mx,my,markermextc] = vis_markers((xlims-ox)*1000,(ylims-oy)*1000,npx,npy,mx,my,100*mextc);
    
    mx = mx + ox*1000;
    my = my + oy*1000;
    
    if j == 1
        xlimits = [mx(1) mx(end)];
        ylimits = [my(1) my(end)];
    end
    
    % Plot the material proporties
    clf
    pcolor(mx/1000,my/1000,markermextc);    hold on
    if ~isempty(tempc)
        [C,h]=contour(G.gridx/1000,G.gridy/1000,G.tk1-273,tempc,'w','linewidth',1);
    end
    
    colormap(cmap);
    shading interp;
    caxis([0 5])
    colorbar
    
    set(gca,'Ydir','reverse')
    title(['Melt Fraction, Step=',num2str(j),' Myr=',num2str(timesum*1e-6/(yr2sec))]);
    xlabel('x, km')
    ylabel('y, km')
    drawnow
end

if ~isempty(tempc)
    clabel(C,h,'manual','rotation',0)
end

end

