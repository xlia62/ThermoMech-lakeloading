function cmap = markercolors(nc)
% Colormap  with labels  

cmap{1} = zeros(nc,3);
% 'Sticky Air';
cmap{1}(1,:) = [255 255 255];
cmap{2}{1} = 'Air';
% 'Sticky Water';
cmap{1}(2,:) = [179 240 255];
cmap{2}{2} = 'Water';
%'Sediments';
cmap{1}(3,:) = [255 153 153];
cmap{2}{3} = 'Sed';
% 'Upper Oceanic Crust';
cmap{1}(4,:) = [102 102 0];
cmap{2}{4} = 'UOC';
% 'Lower Oceanic Crust';
cmap{1}(5,:) = [0 51 0];
cmap{2}{5} = 'LOC';
% 'Upper Continental Crust';
cmap{1}(6,:) = [153, 102, 0];
cmap{2}{6} = 'UCC';
% 'Lower Continental Crust';
cmap{1}(7,:) =  [102 53 0];
cmap{2}{7} = 'LCC';
% 'Oceanic Lithospheric Mantle';
cmap{1}(8,:) =  [51 0 255];
cmap{2}{8} = 'OLM';
% 'Continental Lithospheric Mantle';
cmap{1}(9,:) = [255 204 0];
cmap{2}{9} = 'CLM';
% 'Asthenospheric Mantle';
cmap{1}(10,:) = [255 255 153];
cmap{2}{10} = 'AsthM';
% 'Hydrated Mantle';
cmap{1}(11,:) = [0 115 182];
cmap{2}{11} = 'HydrM';

% Optional rock types
if nc > 11
    % 'Upper Oceanic Plateau';
    cmap{1}(12,:) = [100 100 100];
    cmap{2}{12} = 'UOP';
end

if nc > 12
    % 'Lower Oceanic Plateau';
    cmap{1}(13,:) = [0 0 0];
    cmap{2}{13} = 'LOP';
end

if nc > 13
    % 'Plume Material';
    cmap{1}(14,:) = [120 171 48];%0 102 0];
    cmap{2}{14} = 'Plume';
end

if nc > 14
    % Melt Material
    cmap{1}(15,:) = [255 0 0];
    cmap{2}{15} = 'Melt';
end

if nc > 15
    % Lower mantle
    cmap{1}(16,:) = [173 235 173];
    cmap{2}{16} = 'LM';
end

if nc > 16
    % Eclogite < 50%
    cmap{1}(17,:) = [230 230 230];
    cmap{2}{17} = 'EC';
end

if nc > 17
    % Eclogite > 50%
    cmap{1}(18,:) = [166 166 166];
    cmap{2}{18} = 'EC>50%';
end


cmap{1} = cmap{1}/255;
    