% Define Material Properties, Build Geometry, Construct Boundary Conditions
%
% Edit this file to change:
%
%         material properties
%         model geometry
%         boundary conditions
%

% -------------------------------------------------------------------------
%% Material properties
% -------------------------------------------------------------------------
MouchaMaterials

% -------------------------------------------------------------------------
%% Defining Phase Transformations
% (1) Basalt to eclogite
% (2) 410 km
% (3) 660 km
% -------------------------------------------------------------------------
MPHASE = zeros(8,4);

% Eclogite (prograde and retrograde)

% Possible Phase Transformation types:
% false = no phase transformations
% 'depth' = pressure based according MPHASE(4,:)
% 'thermodynamic' = based on thermodynamic criteria (Gibbs Free Energy)
%                   with blueshicst cut-off
% 'kinetic' = based on kinetics

RxType =  'thermodynamic';
phase_function = 'hyperbollic';   % Univeriate or hyperbollic

MPHASE(1,1:2) = [175e3 175e3];    % Kinetic Activation energy [prograde retrograde]
MPHASE(2,1:2) = [1e-6 1e-6];      % Kinetic Pre-exponential factor [prograde retrograde]
MPHASE(3,1:2) = 4;                % Avrami exponent 
                                  % = 4 means rate of transformation is not constant, but starts 
                                  %     at zero and ends at zero
                                  % = 1 reaction is fast at first than slow
MPHASE(4,1:2)= 1.5e9;             % Fixed Pressure phase transition (Pa)
MPHASE(5,1:2) = 200;             % Density Change (kg/m^3)

% 410 km
MPHASE(4,3) = 1.5;              % Clapeyron slope (MPa/K)
MPHASE(5,3) = 206;              % Density change (kg/m^3)
MPHASE(6,3) = 13.5;             % Reference phase change pressure (GPa)
MPHASE(7,3) = 1810;             % Reference phase change temperature (K)
MPHASE(8,3) = 1;                % Viscosity increase factor
MPHASE(9,3) = 0.1;                % Phase transition range for hyperbolic (0.1 GPa = 20 km)

% 660 km
MPHASE(4,4) = -1.3;             % Clapeyron slope (MPa/K)
MPHASE(5,4) = 322;              % Density change (kg/m^3)
MPHASE(6,4) = 23;               % Reference phase change pressure (GPa)
MPHASE(7,4) = 1940;             % Reference phase change temperature (K)
MPHASE(8,4) = 10;                % Viscosity increase factor
MPHASE(9,4) = 0.1;                % Phase transition range for hyperbolic (0.1 GPa = 20 km)

% -------------------------------------------------------------------------
%% Defining lithological temperature model constants
% -------------------------------------------------------------------------

% Sticky Layer thickness water or air (m)
sticky_layer = 20000;

% Water Level Tracking (m) set to empty if no water
% Water depth = sticky_layer - water_lev
water_lev = 18000;

% continent elevation (m) -- should use isostatic balance
% height above water will be water_lev - cont_lev
cont_lev = 17000;

% Sediment Thickness (m)
sed_thick = 2000;

% Bottom of sediments ocean
bot_osed = sticky_layer + sed_thick;

% Upper Ocean Crust thickness (m)
upper_ocrust_thick = 3000;

% lower Ocean Crust thickness (m)
lower_ocrust_thick = 4000;

% upper continent crust thickness (m)
upper_ccrust_thick = 20000;

% lower continent crust thickness (m)
lower_ccrust_thick = 15000;

% continen mantle lithosphere thick (m)
cont_mantle_lith_thick = 115000;

% Moho depth (m)
moho = (upper_ccrust_thick + lower_ccrust_thick);

% Continental Lithosphere thickness (m)
cont_lith_thick = moho + cont_mantle_lith_thick;

% Not used but needed to set to false
add_plume = false;

%--------------------------------------------------------------------------
%% Initial Temperature parameters
%--------------------------------------------------------------------------

% Moho temperature (K)
tmoho = 600 - tabsolute;

% Thermal diffusivity of the mantle, m^2/s
kappa = 1e-6;    

% Oceanic lithosphere thickness (1100oC isotherm) (m)
ocean_LAB = 1100 - tabsolute;

% Age of Oceanic Plate (s)
ocean_age = 60e6*yr2sec;  % Subducting

% Convert deg-C to K
ttop = ttop - tabsolute;

% Temperature at bottom of model space
mantle_thick = ysize-sticky_layer;
tbottom = tpotential - tabsolute + tgrad*mantle_thick;

% Continent Geotherm piece-wise linear function
%--------------------------------------------------------------------------

% Temperature of adiabat at desired depth
Tl = tpotential - tabsolute + tgrad*cont_lith_thick;

% Crust section -- linear: slope and intercept
mcrust = (tmoho-ttop) / moho;
bcrust = tmoho - mcrust*moho;

% Mantle lithosphere -- linear
% Moho to base of lithosphere, intersect with adiabat: slope and intercept
mlith = (Tl - tmoho)/(cont_lith_thick-moho);
blith = Tl - mlith*cont_lith_thick;

% Oceanic Geotherm -- Half-space cooling geotherm
%--------------------------------------------------------------------------

% Mantle adiabat
z = linspace(0,mantle_thick,200);
T_adiabat = tpotential - tabsolute + tgrad*z;

% Scale half-space cooling with mantle adiabat
T_ocean = ttop + (T_adiabat-ttop).*erf(z./(2*sqrt(kappa*ocean_age)));

% Oceanic lithosphere thickness determined by the ocean_LAB isotherm
k = find(T_ocean>=ocean_LAB,1);
if isempty(k)
    error('Model space exceeds isoLAB')
else
    oceanic_lith_thick = z(k);
end

% Lithosphere asthenosphere boundary (m) 
oceanic_lith_base = sticky_layer + oceanic_lith_thick;
cont_lith_base = cont_lev + cont_lith_thick;

% -------------------------------------------------------------------------
%% Defining lithological structure of the model
% -------------------------------------------------------------------------

% Marker counter
mm1=0;
rng('default') % seed

% Initial trench position (m) -- middle of model space
TrenchX = 1900e3;
% Initial dip angle
subd_angle = 30;         % (deg)
% Width of weak zone
weak_zone_width = 10000; % (m)

% Find slope and intercept for the equation of the two bounding lines:
% y1 = mx + b1,  y2 = mx + b2
y0 = sticky_layer;
m =  tand(subd_angle);
b2 = y0 - m*TrenchX;  % Left bound intercept
b1 = b2 + m*weak_zone_width;  % Right bound intercept

% Loop over makers initial loop
for mm1 = 1:marknum
    
    %--------------------------------------------------------
    % Initial rock distribution
    %--------------------------------------------------------
        
    % Asthenosphere
    MI(mm1) = 10;
    
    % Sticky air
    if(MY(mm1) <= sticky_layer)
        MI(mm1)= 1;
    end
    
    % Sticky water
    if ~isempty(water_lev)
        if (MY(mm1) >= water_lev  && MY(mm1) <= sticky_layer)
            MI(mm1) = 2;
        end
    end
    
    % Sediments
    % Bottom of sediments
    bot_sed = sed_thick + sticky_layer;
    if(MY(mm1)>sticky_layer && MY(mm1)<= bot_sed)
        MI(mm1)=3;
    end    
   
    % Basaltic crust
    % Bottom of upper oceanic crust
    bot_ocrust1 = bot_osed + upper_ocrust_thick;
    if(MY(mm1)>bot_osed && MY(mm1)<=bot_ocrust1 )
        MI(mm1)=4;
    end
    
    % Gabbroic crust
    % Bottom of lower oceanic crust
    bot_ocrust2 = bot_ocrust1 + lower_ocrust_thick;
    if(MY(mm1)>bot_ocrust1 && MY(mm1)<=bot_ocrust2)
        MI(mm1)=5;
    end
        
    %--------------------------------------------------------
    % Subducting Old Ocean Plate (left)
    %--------------------------------------------------------
    
    % Mantle lithosphere
    if(MY(mm1)>bot_ocrust2 && MY(mm1)<=oceanic_lith_base)
        MI(mm1) = 8;
    end
    
    %--------------------------------------------------------
    % Weak subduction zone
    %--------------------------------------------------------
    
    % Location relative to subduction zone, e.g if d1 < 0 right of, d1 > 0 left of
    d1 = m*MX(mm1) + b1 - MY(mm1);
    d2 = m*MX(mm1) + b2 - MY(mm1);    
    
    if d1 >= 0 && d2 <= 0 && MY(mm1) >= sticky_layer && MY(mm1) <= oceanic_lith_base
        MI(mm1) = 11;
    end
    
    %--------------------------------------------------------
    % Over-riding Continent Plate (right)
    %--------------------------------------------------------
    
    if d2 >= 0
        
        bot_ccrust1 = cont_lev + upper_ccrust_thick;
        % Upper continental crust starting at continent_level
        if (MY(mm1)>cont_lev && MY(mm1)<=bot_ccrust1)
            MI(mm1) = 6;
        end
        
        bot_ccrust2 = bot_ccrust1 + lower_ccrust_thick;
        % Lower continental crust starting at upper continent depth
        if (MY(mm1)>bot_ccrust1 && MY(mm1)<=bot_ccrust2)
            MI(mm1) = 7;
        end
        
        cont_LAB = bot_ccrust2 + cont_mantle_lith_thick;
        % Continental Mantle lithospre
        if(MY(mm1)>bot_ccrust2 && MY(mm1)<=cont_LAB)
            MI(mm1) = 9;
        end
        
    end
    
end


%--------------------------------------------------------------------------
% Initial temperature structure
%--------------------------------------------------------------------------

for mm1 = 1:marknum   
    
    % Adiabatic Temperature Gradient in the Mantle
    MTK(mm1)=tbottom-tgrad*(ysize-MY(mm1));
    
    % Sticky air or sticky water
    if(MI(mm1) == 1 || MI(mm1) == 2)
        MTK(mm1) = ttop;
    end
    
    % Oceanic geotherm for subducting plate of specified age
    % after Turcotte & Schubert (2002)
    d2 = m*MX(mm1) + b2 - MY(mm1);  % Include weak zone

    if( MY(mm1) > sticky_layer && d2 <= 0 )
       z = MY(mm1) - sticky_layer;             
       MTK(mm1) = ttop + (MTK(mm1)-ttop).*erf(z/(2*sqrt(kappa*ocean_age)));
    end
    
    % Continental geotherm for overriding plate
    % Linear to moho temperature
    if MI(mm1) == 6 || MI(mm1) == 7
        z = MY(mm1) - cont_lev;
        MTK(mm1) = mcrust*z + bcrust;
    end
    
    if MI(mm1) == 9
        z = MY(mm1) - cont_lev;
        MTK(mm1) = mlith*z + blith;
    end
    
end

% Save Number of markers marknum=mm1
marknum = mm1;

% -------------------------------------------------------------------------
%% Initial elevation for topography profile
% -------------------------------------------------------------------------

% Intial elevation at bottom of sticky layer
gridt(2,:) = sticky_layer;

% Set initial continent elevation
k = gridt(1,:) >= TrenchX;
gridt(2,k) = cont_lev;

% Save initial topography
% topotime(1,1) = 0;
% topohigh(1,:) = gridt(2,:);
% if ~isempty(water_lev)
%     topo_lev=zeros(ntimestep+1,tnum);
%     topowater(1,1) = water_lev;
% end