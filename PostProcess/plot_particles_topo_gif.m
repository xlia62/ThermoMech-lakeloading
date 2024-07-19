function plot_particles_gif(input_dir,nframes,xlims,ylims)

%input_dir = 'Africa2019/lake1200/Tmoho600_dT0';
output_movie = input_dir;
output_image = input_dir;
%nframes = 100;
sticky_layer = 20;
%xlims = [150 225];
%ylims = [-1 20];
marker_size = 20;


yr2sec= 365.25*24*3600;

% colors
c1 = [0.6350, 0.0780, 0.1840];
c2 = [0.8500, 0.3250, 0.0980];
c3 = [139 69 19]/255;
flt = [0.4940, 0.1840, 0.5560];
wtr = [0 0.447 0.7410];
sds = [0.9290 0.6940 0.1250];
mtl = [153 204 0]/255;
ast = [.47 .67 .19];
bslt = [47 79 79]/255;

% figure('Position',[66 250 1167 385],'color','white')
%figure('Position',[17 130 1262 350],'color','white')
count = 0;
ncount = nframes;
for i = 1:nframes
    
    clf
    
    count = count + 1;
    load([input_dir,'/markers_',num2str(i),'.mat'])
    load([input_dir,'/grids_',num2str(i),'.mat'])
    
    subplot(211)
    [X,Y]=meshgrid(gridt(1,:)/1000,linspace(0,FastScape.yl,FastScape.ny));
    Z = FastScape.topo+100;
    surf(X,Y,Z)
    demcmap('inc',Z,100)
    colorbar
    view(-8,70)
    shading interp
    lightangle(-45,30)
    h.FaceLighting = 'gouraud';
    h.AmbientStrength = 0.3;
    h.DiffuseStrength = 0.8;
    h.SpecularStrength = 0.9;
    h.SpecularExponent = 25;
    h.BackFaceLighting = 'unlit';
    material dull
    
    subplot(212)
    mx = MX/1000;
    my = MY/1000 - sticky_layer;
    
    k = mx >= xlims(1) & mx <= xlims(2) & my>= ylims(1) & my<= ylims(2);
    mx = mx(k);
    my = my(k);
    mi = MI(k);
    
    crust1 = mi == 4;
    crust2 = mi == 6;
    fault = mi == 11;
    water = mi == 2;
    seds = mi == 3;
    crust3 = mi == 7;
    mantle = mi == 9;
    astheno = mi == 10;
    basalt = mi == 12;  % Basalt

    scatter(mx(astheno),my(astheno),marker_size,'MarkerFaceColor',ast,'MarkerEdgeColor',ast,...
        'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)
    hold on
    scatter(mx(mantle),my(mantle),marker_size,'MarkerFaceColor',mtl,'MarkerEdgeColor',mtl,...
        'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)
    scatter(mx(crust3),my(crust3),marker_size,'MarkerFaceColor',c3,'MarkerEdgeColor',c3,...
        'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)
    scatter(mx(crust2),my(crust2),marker_size,'MarkerFaceColor',c2,'MarkerEdgeColor',c2,...
        'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)
    scatter(mx(crust1),my(crust1),marker_size,'MarkerFaceColor',c1,'MarkerEdgeColor',c1,...
        'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)
    scatter(mx(basalt),my(basalt),marker_size,'MarkerFaceColor',bslt,'MarkerEdgeColor',bslt,...
        'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)
    scatter(mx(fault),my(fault),marker_size,'MarkerFaceColor',flt,'MarkerEdgeColor',flt,...
        'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)
    scatter(mx(seds),my(seds),marker_size,'MarkerFaceColor',sds,'MarkerEdgeColor',sds,...
        'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)
    scatter(mx(water),my(water),marker_size,'MarkerFaceColor',wtr,'MarkerEdgeColor',wtr,...
        'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)
    
%   axis image
%   plot isotherms
    contour(G.gridx/1000,G.gridy/1000 - sticky_layer,G.tk1-273,[50 150],'linewidth',1,'color','k')
    ylim(ylims)
    xlim(xlims)
    set(gca,'Ydir','reverse')
    set(gca,'FontSize',14)
    
%   Plot horizons
    modeltime = timesum/yr2sec/1e6;
    nh = floor(modeltime)+1;
    for i = 1:nh
        plot(gridt(1,:)/1000,H(i,:)/1000-sticky_layer,'linewidth',1,'color','b')
    end
    plot(gridt(1,:)/1000,gridt(2,:)/1000-sticky_layer,'linewidth',1,'color','k')

    title([num2str(modeltime),' My'],'HorizontalAlignment','right','Position',[200 -1 0])

    f = getframe(gcf);
    im = frame2im(f);
    
    [A,map] = rgb2ind(im,256,'nodither');
    if count == 1
        imwrite(A,map,[output_movie,'.gif'],'gif','LoopCount',Inf,'DelayTime',.9);
    else
        imwrite(A,map,[output_movie,'.gif'],'gif','WriteMode','append','DelayTime',.9);
    end
%     
%     if count == 0
%         [im,map] = rgb2ind(f.cdata,256,'nodither');
%         im(1,1,1,ncount) = 0;
%     else
%     im(:,:,1,count+1) = rgb2ind(f.cdata,map,'nodither');   
%     end
end

% Add small delay at end before restart
imwrite(A,map,[output_movie,'.gif'],'gif','WriteMode','append','DelayTime',.1);

%print('-dtiff','-r200',[output_image,'.tiff'])
