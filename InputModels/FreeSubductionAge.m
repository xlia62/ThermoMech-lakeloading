% Define Material Properties, Build Geometry, Construct Boundary Conditions
%

%% Define Material Properties by loading a file

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
MPHASE(5,1:2) = 350;             % Density Change (kg/m^3)

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
MPHASE(8,4) = 1;                % Viscosity increase factor
MPHASE(9,4) = 0.1;                % Phase transition range for hyperbolic (0.1 GPa = 20 km)

% -------------------------------------------------------------------------
%% Defining lithological temperature model constants
% -------------------------------------------------------------------------

% Sticky Layer thickness water or air (m)
sticky_layer = 20000;

% Water Level Tracking (m) set to empty if no water
% Water depth = sticky_layer - water_lev
water_lev = 14000;

% Sediment Thickness (m)
sed_thick = 2000;

% Bottom of sediments ocean
bot_osed = sticky_layer + sed_thick;

% Upper Ocean Crust thickness (m)
upper_ocrust_thick = 3000;

% lower Ocean Crust thickness (m)
lower_ocrust_thick = 4000;

% Total oceanic crust thickness (m)
ocrust_thick = upper_ocrust_thick + lower_ocrust_thick;

% Weak zonde thickness (m)
wzone_thick = 5e3;

% Trench position (m)  for 5000km use 3000, for 3000 use 1800
trench_x = 1800e3;

% Age of Oceanic Plate (Ma)
plate_age_sp = 60;  % Subducting

% Age of Oceanic Plate (Ma)
plate_age_op = 40; % Overriding

% Length of subducting plate (km)
length_plate_sp = trench_x/1000;

% Length of overriding plate (km)
length_plate_op = (xsize-trench_x)/1000;

% Add plume at given time (Ma)
add_plume = false;

%--------------------------------------------------------------------------
%% Initial Temperature parameters
%--------------------------------------------------------------------------

% Thermal diffusivity of the mantle, m^2/s
kappa = 1e-6;    

% Oceanic lithosphere thickness (1000oC isotherm) (m)
isoLAB = 1100-tabsolute;

% Bottom of thermal lithosphere
ThermLith = 1300-tabsolute;

% Convert deg-C to K
ttop = ttop - tabsolute;

% Temperature at bottom of model space
mantle_thick = ysize-sticky_layer;
tbottom = tpotential - tabsolute + tgrad*mantle_thick;

% Thermal lithosphere
ThermLAB = 1100 - tabsolute;

% Compute subduction temperatures
theta = 1:55;
rcurve = 250e3;  % radius of curvature in m

% Find all markers that will be used to for oceanic plate thermal structure 
mm1 = MY > sticky_layer & (MY-sticky_layer) <= rcurve;

% Assumes the following, a 60 Ma subductiong plate and a 30 Ma overriding
% plate --- needs to be made more general.
MTK(mm1) = OceanOceanSubThermal(MX(mm1)/1000,(MY(mm1)-sticky_layer)/1000,plate_age_sp,plate_age_op,length_plate_sp,length_plate_op,theta,tpotential,rcurve/1000,tgrad) - tabsolute;

% -------------------------------------------------------------------------
%% Defining lithological structure of the model
% -------------------------------------------------------------------------

% Crust zone polygon
crust_radius = rcurve - ocrust_thick;

x1 = trench_x + crust_radius*cosd(90-theta-3);
x2 = trench_x + rcurve*cosd(90-theta-3);
y1 = crust_radius - crust_radius*sind(90-theta-3) + sticky_layer;
y2 = crust_radius - rcurve*sind(90-theta-3) + sticky_layer;
k = y1 < sticky_layer+sed_thick;
x1(k) = [];
y1(k) = [];
k = y2 < sticky_layer+sed_thick;
x2(k) = [];
y2(k) = [];

xpc = [x1 fliplr(x2)];
ypc = [y1 fliplr(y2)];
   
% Weak zone polygon
wzone_radius = rcurve;

x1 = trench_x + wzone_radius*cosd(90-theta-3);
x2 = trench_x + (wzone_radius+wzone_thick)*cosd(90-theta-3);
y1 = wzone_radius - wzone_radius*sind(90-theta-3) + sticky_layer - ocrust_thick;
y2 = wzone_radius - (wzone_radius+wzone_thick)*sind(90-theta-3) + sticky_layer - ocrust_thick;
k = y1 < sticky_layer+sed_thick+ocrust_thick;
x1(k) = [];
y1(k) = [];
k = y2 < sticky_layer+sed_thick+ocrust_thick;
x2(k) = [];
y2(k) = [];

xpw = [x1 fliplr(x2)];
ypw = [y1 fliplr(y2)];

% Loop over makers to set material properties and intial thermal structure
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
    
    %--------------------------------------------------------
    % Thermal models
    %--------------------------------------------------------
    
    % Adiabatic Temperature Gradient in the Mantle
    if (MY(mm1)-sticky_layer) > rcurve
        MTK(mm1)=tbottom-tgrad*(ysize-MY(mm1));
    end
    
    % Sticky air or sticky water
    if(MI(mm1) == 1 || MI(mm1) == 2)
        MTK(mm1) = ttop;
    end
    
    %--------------------------------------------------------
    % Oceanic Plates defined by temperature
    %--------------------------------------------------------
    
    % Mantle lithosphere
    if(MY(mm1)>sticky_layer && MTK(mm1) <= ThermLAB)
        MI(mm1) = 8;
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
end

% Create crust zone
[in,on]=inpolygon(MX,MY,xpc,ypc);
MI(in|on) = 5;

% Create weak zone
[in,on]=inpolygon(MX,MY,xpw,ypw);
MI(in|on) = 11;


% -------------------------------------------------------------------------
%% Initial elevation for topography profile
% -------------------------------------------------------------------------

% Intial elevation at bottom of sticky layer
gridt(2,:) = sticky_layer;

% Set initial continent elevation
%k = gridt(1,:) >= ix & gridt(1,:);
%gridt(2,k) = cont_lev;

% Save initial topography
% topotime(1,1) = 0;
% topohigh(1,:) = gridt(2,:);
% if ~isempty(water_lev)
%     topo_lev=zeros(ntimestep+1,tnum);
%     topowater(1,1) = water_lev;
% end