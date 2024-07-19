function animate_contour(inpath,nbegin, nfreq,nend,xlims,ylims, MPROP)

% Inpath = folder where grid_###.mat and markers_###.mat files exist
% nfreq  = frequency of time stpes when files are saved
% nbegin = start number of time steps 
% nend = end number of time steps
% npx,npy = number of image pixels in x and y direction
% cmap = colormap
% tempc = temperature cont

%ours, ommit adding contours if empty
npx = 400;
npy = 200;
cmap = jet;
tempc = 5;


k = [nbegin:nfreq:nend];

yr2sec = 365.25*24*3600;

figure
pos = get(gcf,'Position');
width = pos(3) * 3;
height = pos(4) * 1.5;
set(gcf,'Position',[pos(1) pos(2) width height]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Optional: write to video


figname = split(inpath, '/');
video_name = ['output/AfricaModels2022/Figure/extraction/',char(figname(end)),'_',num2str(nbegin),'-', num2str(nend), '_', MPROP];
%video_name = 'test';
video = VideoWriter(video_name);
try
    close(video);
catch
    'here'
end
video.FrameRate = 5;   % if every timestep is saved, a little fast
open(video);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%myfig = figure('Position',[17 130 1262 489],'color','white')


for j = k
    
    % Load files
    infile = [inpath,'/markers_',num2str(j),'.mat'];
    load(infile);
    infile = [inpath,'/grids_',num2str(j),'.mat'];
    load(infile);
    
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
            unit = 'Cummulative Melt Fraction';
        case 'MEXTC' % Cummulative Extracted Melt Fraction
            MGII = MEXTC;
            unit = 'Cummulative Extracted Melt Fraction';
        case 'MSXY' % Cummulative Extracted Melt Fraction
            MGII = MSXY;
            unit = 'shear stress (Pa)';
        case 'MSXX' % Cummulative Extracted Melt Fraction
            MGII = MSXX;
            unit = 'normal stress (Pa)';
    end
    
    % Replace Air or Water with NaN
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
    
    
    [mx,my,markerrho] = vis_markers((xlims-ox)*1000,(ylims-oy)*1000,npx,npy,mx,my,mgii);
    
    mx = mx + ox*1000;
    my = my + oy*1000;
    
    if j == 1
        xlimits = [mx(1) mx(end)];
        ylimits = [my(1) my(end)];
    end
    
    % Plot the material proporties
    clf;
    pcolor(mx/1000,my/1000,markerrho);
    colors1 = [.7 .7 .7; .4 .4 .4; .005 .005 .005];
    n = 50;
    m = size(colors1,1);
    t0 = linspace(0,1,m)';
    t = linspace(0,1,n)';
    r = interp1(t0,colors1(:,1),t);
    g = interp1(t0,colors1(:,2),t);
    b = interp1(t0,colors1(:,3),t);
    cmap2 = [r,g,b];
    colormap(jet);
    shading interp;
%    set(gca,'renderd
    
    set(gca,'Ydir','reverse')
    title(['Step=',num2str(j),' Myr=',num2str(timesum*1e-6/(yr2sec))],'FontSize', 24);
    xlabel('Horizontal Position (km)','FontSize', 20,'FontWeight','bold')
    ylabel('Depth (km)', 'FontSize',20,'FontWeight','bold')
    c = colorbar;
    ylabel(c,unit,'FontSize',20);
    axis equal
    drawnow ();
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % For video
    %pause(.02)
    f = getframe(gcf);
    video.writeVideo(f);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for video
close(video)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

