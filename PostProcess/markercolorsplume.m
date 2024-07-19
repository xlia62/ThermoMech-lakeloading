function cmap = markercolorsplume(nc)
% Colormap
cmap = zeros(nc,3);
% 'Sticky Air';
cmap(1,:) = [255 255 255];
% 'Sticky Water';
cmap(2,:) = [0 204 255];
%'Sediments';
cmap(3,:) = [255 94 246];
% 'Upper Oceanic Crust';
cmap(4,:) = [102 102 0];
% 'Lower Oceanic Crust';
cmap(5,:) = [0 51 0];   
% 'Upper Continental Crust';
cmap(6,:) = [255 153 102];
% 'Lower Continental Crust';
cmap(7,:) =  [102 53 0];
% 'Oceanic Lithospheric Mantle';
cmap(8,:) =  [51 0 255];
% 'Continental Lithospheric Mantle';
cmap(9,:) = [255 204 0];
% 'Asthenospheric Mantle';
cmap(10,:) = [255 255 153];
% 'Hydrated Mantle';
cmap(11,:) = [102 255 102];
if nc > 11
% 'Upper Oceanic Plateau';
cmap(12,:) = [178 178 178];
% 'Lower Oceanic Plateau';
cmap(13,:) = [0 0 0];
end
if nc > 13
% 'Plume Material';
cmap(14,:) = [0 102 0];
end
cmap = cmap/255;


    