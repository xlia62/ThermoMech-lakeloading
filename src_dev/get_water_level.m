function water_level = get_water_level(gridt,water_depth,sticky_layer)

% Uses the minimum point in the topography to set water level according to
% specified water depth

k = islocalmax(gridt(2,:));  
[maxy,iy] = max(gridt(2,k)); % lowest point, since y is depth
k = find(k);
iy = k(iy);

k = find(islocalmin(gridt(2,:)));  % all high points, peaks
ix1 = find(k<iy,1,'last');
ix2 = find(k>iy,1,'first');

% miny = max([gridt(2,k(ix1)) gridt(2,k(ix2)) sticky_layer]);
% 
% depth = min(maxy - miny,water_depth);

water_level = max(maxy-water_depth,sticky_layer);

