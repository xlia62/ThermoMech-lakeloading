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
% GarelMaterials
% GeryaMaterials
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
MPHASE(8,4) = 10;                % Viscosity increase factor
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

% Thermal diffusivity of the mantle, m^2/s
kappa = 1e-6;    

% Oceanic lithosphere thickness (1000oC isotherm) (m)
isoLAB = 1100-tabsolute;

% Bottom of thermal lithosphere
ThermLith = 1300-tabsolute;

% Age of Subdicting Oceanic Plate (s)
plate_age_sp = 60e6*yr2sec;  % Subducting

% Age of Overriding Oceanic Plate (s)
plate_age_op = 40e6*yr2sec; % Overriding

% Add plume at specified time during simulation (Ma) 
add_plume = false;

%--------------------------------------------------------------------------
%% Initial Temperature parameters
%--------------------------------------------------------------------------

% Convert deg-C to K
ttop = ttop - tabsolute;

% Temperature at bottom of model space
mantle_thick = ysize-sticky_layer;
tbottom = tpotential - tabsolute + tgrad*mantle_thick;

% Half-space cooling geotherm
z = linspace(0,mantle_thick,200);
T_adiabat = tpotential - tabsolute + tgrad*z;
T_sp = ttop + (T_adiabat-ttop).*erf(z./(2*sqrt(kappa*plate_age_sp)));
T_op = ttop + (T_adiabat-ttop).*erf(z./(2*sqrt(kappa*plate_age_op)));

% Seismic oceanic lithosphere thickness at isoLAB
k = find(T_sp>=isoLAB,1);
if isempty(k)
    error('Model space exceeds isoLAB')
else
    oceanic_lith_thick_sp = z(k);
end
k = find(T_op>=isoLAB,1);
if isempty(k)
    error('Model space exceeds isoLAB')
else
    oceanic_lith_thick_op = z(k);
end

% Thermal oceanic lithosphere depth at ThermLith
k = find(T_sp>ThermLith,1);
if isempty(k)
    error('Model space exceeds ThermLith')
else
    oceanic_lith_base = z(k);
end

% Bottom of sediments ocean
bot_osed = sticky_layer + sed_thick;

% Upper Ocean Crust thickness (m)
upper_ocrust_thick = 3000;

% lower Ocean Crust thickness (m)
lower_ocrust_thick = 4000;

% Ocean lithosphere asthenosphere boundary (m) 
oceanLAB_sp = sticky_layer + oceanic_lith_thick_sp;
oceanLAB_op = sticky_layer + oceanic_lith_thick_op;
oceanic_lith_base = sticky_layer + oceanic_lith_base;

% -------------------------------------------------------------------------
%% Defining lithological structure of the model
% -------------------------------------------------------------------------

% Marker counter
mm1=0;
rng('default') % seed

% Initial trench position (m) -- middle of model space
TrenchX = 1500e3;

% Dip-angle and thickness of weekzone (deg,m)
weakz_thick = 20e3;
weakz_dip = 30;

% Equation describing the weak zone x-position as function of y
weakz_slope = tand(weakz_dip);
weakz_intercept = sticky_layer - weakz_slope*TrenchX;
weakz_x = @(y) (y-weakz_intercept)./weakz_slope;

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
    % Subducting "Old" Ocean Plate (left)
    %--------------------------------------------------------
    
    % Mantle lithosphere
    if(MY(mm1)>bot_ocrust2 && MY(mm1)<=oceanLAB_sp && MX(mm1)<=weakz_x(MY(mm1)))
        MI(mm1) = 8;
    end
    
    %--------------------------------------------------------
    % Overriding "Young" Ocean Plate
    %--------------------------------------------------------
        
    % Mantle lithospre
    if(MY(mm1)>bot_ocrust2 && MY(mm1)<=oceanLAB_op && MX(mm1)>=weakz_x(MY(mm1)))
        MI(mm1) = 8;
    end
    
    % Weak zones
    % Subduction channel
    if (MY(mm1)> bot_ocrust2 && MY(mm1)<=oceanLAB_sp && MX(mm1)>=weakz_x(MY(mm1)) && MX(mm1) <= weakz_x(MY(mm1)) + weakz_thick)
        MI(mm1) = 11;
    end
    
%     % Subducting plate's ridge
%     if (MY(mm1)> bot_ocrust2 && MY(mm1)<=oceanLAB_sp && MX(mm1)<=oceanic_lith_thick_sp)
%         MI(mm1) = 11;
%     end
%     
%     % Overriding plate's ridge
%     if (MY(mm1)> bot_ocrust2 && MY(mm1)<=oceanLAB_op && MX(mm1)>=xsize0 - oceanic_lith_thick_op)
%         MI(mm1) = 11;
%     end
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
    if( MY(mm1) > sticky_layer && MY(mm1)<=oceanic_lith_base && MX(mm1) <= weakz_x(MY(mm1)))
       z = MY(mm1) - sticky_layer;
       MTK(mm1) = ttop + (MTK(mm1)-ttop).*erf(z/(2*sqrt(kappa*plate_age_sp)));
    end
    
    % Oceanic geotherm for overriding plate of specified age and weak zone
    % after Turcotte & Schubert (2002)
    if( MY(mm1) > sticky_layer && MY(mm1)<=oceanic_lith_base && MX(mm1) > weakz_x(MY(mm1)))   
        z = MY(mm1) - sticky_layer;
        MTK(mm1) = ttop + (MTK(mm1)-ttop).*erf(z/(2*sqrt(kappa*plate_age_op)));
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
%k = gridt(1,:) >= ix & gridt(1,:);
%gridt(2,k) = cont_lev;

% Save initial topography
% topotime(1,1) = 0;
% topohigh(1,:) = gridt(2,:);
% if ~isempty(water_lev)
%     topo_lev=zeros(ntimestep+1,tnum);
%     topowater(1,1) = water_lev;
% end