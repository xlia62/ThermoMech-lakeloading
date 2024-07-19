function [etime,vext,min_lith,min_crust]=thinning_factors(inpath,nfreq,ntotal,xlims,lith_depth,continental)

% Estimate and plot thinning factors for crust and lithosphere

if nfreq > 1
    k = [1 nfreq:nfreq:ntotal];
else
    k = [nfreq:nfreq:ntotal];
end

nframes = length(k);

%k([49 53 57 58]) = [];

yr2sec = 365.25*24*3600;
etime = zeros(nframes,1);
vext = etime;

% Crust markers, continental == 1
upper_crust = 6;
lower_crust = 7;

% Oceanic crust
if ~continental
    upper_crust = 4;
    lower_crust = 5;
end

% Generate bins of 2 km
edges = xlims(1):2:xlims(2);
nbins = length(edges)-1;

xlims = xlims * 1000;
lith_depth = lith_depth*1000;

for j = k
    
    load([inpath,'/markers_',num2str(j),'.mat'],'MI','MX','MY')
    load([inpath,'/grids_',num2str(j),'.mat'])
    
    etime(j) = timesum/yr2sec/1e6;       % Myr
    vext(j) = (G.vx1(1,end)-G.vx1(1,1))*yr2sec*1000;  % mm/yr
    
    % Restrict to 400-600 km (100 km on either side of center)
    % Crust (upper and lower)
    c1 = (MI == upper_crust | MI == lower_crust | MI == 11) & MX >= xlims(1) & MX <= xlims(2);
    c1x = MX(c1)/1000;
    c1y = MY(c1)/1000;
    
    % Lithosphere (top of asthenosphere)
    c2 = MI == 10 & MX >= xlims(1) & MX <= xlims(2) & MY < lith_depth;
    c2x = MX(c2)/1000;
    c2y = MY(c2)/1000;
    
    i1 = discretize(c1x,edges);
    i2 = discretize(c2x,edges);
    
    t1 = zeros(nbins,1);
    t2 = t1;
    t3 = t1;
    
    for i = 1:nbins
        k1 = i1 == i;
        t1(i) = max(c1y(k1))-min(c1y(k1));
        k2 = i2 == i;
        t2(i) = min(c2y(k2))- min(c1y(k1));
    end
    
    % Thinnest crust and lithosphere
    min_lith(j) = min(t2);
    min_crust(j) = min(t1);
    
end

