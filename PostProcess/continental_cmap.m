function cmap=continental_cmap(opt)

% Air
cmap =[255 255 255];

% Sediments255 94 246;
% New Crust
% Upper Crust 
% Lower Crust
% Lithosphere
% Asthenosphere
cmap = [cmap;
    102 102 0;
    0 51 0;   
    255 153 102;
    102 53 0;
    102 102 102;
    102 255 102;
    204 0     0;
    255 204   0];

% Add plume material
if opt
    cmap = [cmap; 
        0 102 51;
        102 0 102];
end

cmap = cmap./255;