function fig = animate_matprop_remote(inpath,nfreq,ntotal,npx,npy,tempc,xlims,ylims)

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
    fig(j).file = infile;
    fig(j).time = timesum*1e-6/yr2sec;
    fig(j).mx = mx/1000;
    fig(j).my = my/1000;
    fig(j).markcom = markcom;
    
    if ~isempty(tempc)
        fig(j).tempc = tempc;
        fig(j).gridx = G.gridx/1000;
        fig(j).gridy = G.gridy/1000;
        fig(j).tk1 = G.tk1-273;
    end
    
end

end

