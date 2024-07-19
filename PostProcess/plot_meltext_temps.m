function D=plot_meltext_temps(inpath,nfreq,ntotal,xlims,ylims)

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

c = 0;
D = [];

for j = k
    
    % Load files
    infile = [inpath,'/markers_',num2str(j),'.mat'];
    chk = whos('-file',infile);
    
    if ismember('MEXT',{chk.name})
        load(infile,'MX','MY','MI','MTK','MEXT','MPR');
        no_meltext = false;
    else
        no_meltext = true;
    end
    
    infile = [inpath,'/grids_',num2str(j),'.mat'];
    load(infile);
    
    if ~isempty(xlims)
       i = MX/1000 >= xlims(1) & MX/1000 <= xlims(2);
       mi = MI(i);
       mx = MX(i);
       my = MY(i);
       mtemp = MTK(i);
       mpr = MPR(i);
       mext = MEXT(i);
    else
       xlims = [G.gridx(1) G.gridx(end)]/1000;
       mi = MI
       mx = MX;
       my = MY;
       mtemp = MTK;
       mpr = MPR;
       mext = MEXT;

    end
    if ~isempty(ylims)
        i = my/1000 >= ylims(1) & my/1000 <= ylims(2);
        mi = mi(i);
        my = my(i);
        mx = mx(i);
        mtemp = mtemp(i);
        mpr = mpr(i);
        mext = mext(i);
    else
        ylims = [G.gridy(1) G.gridy(end)]/1000;
    end
    
    % plot only extacted melt temperatures
    i = mext > 0;
    
    if sum(i) > 5
        c = c + 1;
        
        D(c).mtime = timesum*1e-6/(yr2sec);   
        D(c).lims = [xlims ylims];
        D(c).mx = mx(i);
        D(c).my = my(i);
        D(c).mi = mi(i);
        D(c).mtemp = mtemp(i);
        D(c).mpr = mpr(i);
        D(c).mext = mext(i);
    end
    
end

end

