function animate_emplacement(inpath,nfreq,ntotal,xlims,ylims,res,time_units,sticky_layer)


if nfreq > 1
    k = [1 nfreq:nfreq:ntotal];
else
    k = [nfreq:nfreq:ntotal];
end

yr2sec = 365.25*24*3600;

count = 0;
ncount = length(k);

switch time_units
    
    case 'kyr'
        
        yr2sec = yr2sec*1e3;
        time_units = ' kyr';
        
    case 'Myr'
        
        yr2sec = yr2sec*1e6;
        time_units = ' Myr';
        
    otherwise
        
        % do nothing in yrs
        time_units = ' yr';
end

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
    
    infile = [inpath,'/grids_',num2str(j),'.mat'];
    load(infile);
    
    MY = MY -sticky_layer;
    
    % Crop to desired area
    if ~isempty(xlims)
        i = MX/1000 >= xlims(1) & MX/1000 <= xlims(2);
        mx = MX(i);
        my = MY(i);
        mi = MI(i);
        mxm = MXM(i);
    else
        xlims = [G.gridx(1) G.gridx(end)]/1000;
        ox = 0;
        mx = MX;
        my = MY;
        mi = MI;
        mxm = MXM;
    end
    
    if ~isempty(ylims)
        i = my/1000 >= ylims(1) & my/1000 <= ylims(2);
        my = my(i);
        mx = mx(i);
        mi = mi(i);
        mxm = mxm(i);
    else
        ylims = [G.gridy(1) G.gridy(end)]/1000;
    end
    
    setup_patches(mx,my,mxm,mi,res)
    set(gca,'Ydir','reverse','FontSize',12,'color',[0.9 0.9 0.9])
    set(gcf,'color',[0.9 0.9 0.9])
    title(['Material Properties, Time:',num2str(round(timesum/(yr2sec),2)),time_units],'FontSize', 16);
    xlabel('Horizontal Position (km)','FontSize', 12,'FontWeight','bold')
    ylabel('Depth (km)', 'FontSize',12,'FontWeight','bold') 
    set(gca,'Xlim',xlims)
    set(gca,'Ylim',ylims)
    
end

end

function setup_patches(MX,MY,MXM,MI,res)

% colors
c1 = [1.0, 0.750, 0.50];
c2 = [0.90, 0.45 0.];
c3 = [191 128 64]/255;
flt = [0.4940, 0.1840, 0.5560];
wtr = [0 0.447 0.7410];
sds = [148 148 184]/255;
mtl = [153 204 0]/255;
ast = [.47 .67 .19];
mantle_melt = [230 0 0]/255;
crust_melt = [255, 235, 230]/255;
mantle_solid = [77 0 0]/255;



% upper-crust 1
k = MI == 4;

mx = MX(k)'/1000;
if ~isempty(mx)
    my = MY(k)'/1000;
    makepatch(mx,my,res,c1);
end

hold on

% upper-crust 2
k = MI == 6;

mx = MX(k)'/1000;
if ~isempty(mx)    
    my = MY(k)'/1000;
    makepatch(mx,my,res,c2);
end

% sediments
k = MI == 3;

mx = MX(k)'/1000;
if ~isempty(mx)
    my = MY(k)'/1000;
    makepatch(mx,my,res,sds);
end

% lower-crust 1
k = MI == 7;

mx = MX(k)'/1000;
if ~isempty(mx)
    my = MY(k)'/1000;
    makepatch(mx,my,res,c3);
end

% lower-crust 2

% lithosphere
k = MI == 9;

mx = MX(k)'/1000;
if ~isempty(mx)
    my = MY(k)'/1000;
    makepatch(mx,my,2*res,mtl);
end

% asthenosphere
k = MI == 10;

mx = MX(k)'/1000;
if ~isempty(mx)
    my = MY(k)'/1000;
    makepatch(mx,my,2*res,ast);
end

% hydrated mantle melt
k = MI == 11 & MXM > 0;
mx = MX(k)'/1000;
if ~isempty(mx)
    my = MY(k)'/1000;
    makepatch(mx,my,3*res,mantle_melt);
end

% hydrated mantle solid
k = MI == 11 & MXM <= 0;
mx = MX(k)'/1000;
if ~isempty(mx)
    my = MY(k)'/1000;
    makepatch(mx,my,res/2,mantle_solid);
end

% Molten crust
k = MXM > 0 & (MI == 4 | MI == 6);
mx = MX(k)'/1000;
if ~isempty(mx)
    my = MY(k)'/1000;
    makepatch(mx,my,res/2,crust_melt);
end

set(gca, 'Ydir', 'reverse')
axis equal
axis tight
hold off
drawnow

end

function h = makepatch(mx,my,res,c)

X = [mx-res; mx+res; mx+res; mx-res];
Y = [my-res; my-res; my+res; my+res];
h=patch('Xdata',X,'Ydata',Y,'facevertexcdata',c,'facecolor','flat','edgecolor','none','facealpha',0.5);

end

