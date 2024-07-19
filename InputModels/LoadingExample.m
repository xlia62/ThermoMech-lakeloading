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

AfricaMaterials
% The materials have 2 different markers for upper continental crust to
% visualize deformation

% Also there is weak upper continental crust to seed the model with faults

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

RxType =  false; %'thermodynamic';
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
%% Defining lithological model constants
% -------------------------------------------------------------------------

% Sticky Layer thickness water or air (m)
sticky_layer = 20000;

% Water Level Tracking (m) set to empty if no water
% Water depth = sticky_layer - water_lev
water_lev = [];%10750;%20200;

% Sediment Thickness (m)
sed_thick = 0;

% Continent level (m), use this to set continent above or at sea level
cont_lev = sticky_layer;

% Upper Continental Crust thickness (m)
upper_cont_thick = 15000 + sed_thick;  % Thickness includes sediment layer

% How many layers in upper crust to alternate for visualizing deformation
num_crust = 16;

% Lower Continental Crust thickness (m)
lower_cont_thick = 15000;

% Crustal thickness (with sediments) (m)
crust_thick = upper_cont_thick + lower_cont_thick;

% Mantle Lithosphere (m) 
mantle_lith = 105000;

% Lithosphere (LAB) (m) % Specifically, the Thermal lithosphere depth
LAB_depth = mantle_lith + crust_thick;

% Weak Zones

% Size and position of half-disk weak zone (m)
%rweak = 5000;
%xweak = xsize/2;
%yweak = upper_cont_thick  + sticky_layer; % brittle region of mantle lithosphere (model specific)

% Fault
dip_angle = 60;    %(deg)
fault_width = 500; %(m)
x0 = xsize/2;
y0 = sticky_layer;
m = -tand(dip_angle);
b1 = y0 - m*x0;  % Left bound intercept
b2 = b1 + m*fault_width;  % Right bound intercept

% Define load
lx1 = xsize/2 - 10000;
lx2 = xsize/2 + 10000;
lh = 4000;

%--------------------------------------------------------------------------
%% Initial Temperature parameters
%--------------------------------------------------------------------------

% Weak zone temperature perturbation [k];
Tweak = 30;

% Plume Excess Temperature [K]
add_plume = false;  % Ma,  set to false

% Thermal diffusivity of the mantle, m^2/s
kappa = 1e-6;    

% Increase default mantle potential temperature oC
dtpotential = 0;%100;

% Convert deg-C to K
ttop = ttop - tabsolute;
tpotential = tpotential + dtpotential - tabsolute;

% Temperature at bottom of model space
tbottom = 1200 - tabsolute;  % Will be updated below after the

%--------------------------------------------------------------------------
%% Generate initial geotherm and continental strength envelopes
%--------------------------------------------------------------------------

% Continental rifting only -- needs to be modified for general case
% 
% Assemble The geotherm (Temperature Profile to LAB)
zmax = ysize-sticky_layer;    % Depth of geotherm for strength envelop and lithosphere
SE.z = linspace(0,zmax,500)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assemble the geotherm
% n= 0 for piece-wise linear, 1 for exponential heat production 2 for constant heat production
n = 0;

Zl = LAB_depth;    % Mechanical Lithosphere Depth (km)
Ti = ttop + tabsolute; % Ti - Surface Temperature (oC)
Tl = 1200;             % Tl = Temperature of the base of Mantle Lithosphere (oC)

switch n
    
    case 0
        
        % Moho temperature 700oc
        S0 = 600;          % Moho Temperature (oC)
        hr = crust_thick;  % Moho depth (km)
        
        k = [];
        qm = [];
        Zrad = [];

    case 1
        
        k = 2.6;     % Thermal Conductivity
        S0 = 3.0e-6; % Surface heat production rate
        hr = 10000;  % Characteristic drop off - the heat production at z = hr is (1/e) of surface heat production.
        qm = [];
        Zrad = [];

    case 2
        
        k = 2.6;     % Thermal Conductivity
        S0 = 3.0e-6; % Surface heat production rate
        Zrad = 10000; % Constant Radioactive zone% Use for Constant heat production rate
        qm = 0.03; %Heat flow of the continent - ??
        hr = [];
end

% The continetal geotherm
SE.T = ContinentalGeothermForStrengthProfile(n,SE.z,k,hr,S0, qm, Zrad, Ti, Tl, Zl,tgrad,tpotential+tabsolute);
SE.T = SE.T - tabsolute;  % Kelvin

% Calulate the correct bottom temperature
tbottom = SE.T(end);

% Update the thermal gradient to ensure

% Strength Envelopes
SE.eii = 1.6e-15;
% Continental Lithosphere
SE = strengthprofile(sed_thick,MRHO(3,:),MFLOW(3,:),MPL(3,:),...
    upper_cont_thick-sed_thick,MRHO(6,:),MFLOW(6,:),MPL(6,:),...
    lower_cont_thick,MRHO(7,:),MFLOW(7,:),MPL(7,:),...
    mantle_lith,MRHO(9,:),MFLOW(9,:),MPL(9,:),...
    zmax,MRHO(10,:),MFLOW(10,:),MPL(10,:),SE);

save([outpath,'/litho_strength.mat'],'SE'); 

% -------------------------------------------------------------------------
%% Defining lithological structure of the model
% -------------------------------------------------------------------------

% LAB depth with sticky_layer is defined as the depth of 1200oC isotherm
% k = find(SE.T < 1200-tabsolute,1,'last');
% LAB_depth = SE.z(k) + sticky_layer;
LAB_depth = LAB_depth + sticky_layer;

% Marker counter
rng('default') % seed


% Loop over makers
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
    bot_sed = sticky_layer + sed_thick;
    if(MY(mm1)>sticky_layer && MY(mm1)<= bot_sed)
        MI(mm1)=3;
    end
    
    
    % Continental crust
    % Bottom of upper crust
    bot_crust1 = bot_sed + upper_cont_thick;
    if(MY(mm1)>bot_sed && MY(mm1)<=bot_crust1)
        
        MI(mm1) = 6;
        % Layered structure in the sand
        m2 = double(int16(MY(mm1)/(upper_cont_thick/num_crust)-0.5)) + 1;
        m2 = m2 - double(int16(m2/2-0.5))*2;
        if(m2==0)
            MI(mm1)=4;
        end
    end
    
    % Continental crust
    % Bottom of lower crust
    bot_crust2 = bot_crust1 + lower_cont_thick;
    if(MY(mm1)>bot_crust1 && MY(mm1)<=bot_crust2)
        MI(mm1)=7;
    end
    
    % Strong Continetal Mantle Lithosphere
    if(MY(mm1)>bot_crust2 && MY(mm1)<= LAB_depth)
        MI(mm1) = 9;
    end
   
    if(MY(mm1) > 660e3+sticky_layer)
        MI(mm1) = 16;
    end
    
    %-----------
    % load
    %-----------
    
    if MX(mm1) >= lx1 && MX(mm1) <= lx2
        if MY(mm1) >= (sticky_layer - lh) && MY(mm1) <= sticky_layer
            MI(mm1)=4;
        end
    end
           
    %--------------------------------------------------------
    %  Weak zones
    %--------------------------------------------------------
    
    % Thermal perturbation in the base of the lithosphere
    % (to start extension in the middle of the model)
    %     dx = MX(mm1) - xweak;
    %     dy = MY(mm1) - yweak;
    %     dr = (dx^2+dy^2)^0.5;
    %     if (MY(mm1) <= yweak && dr <= rweak)
    %         MI(mm1) = 11;
    % %         MTK(mm1) = MTK(mm1) + Tweak;%*(1-dr/rweak);
    %     end
    
    % Fault zone
    if MI(mm1) == 6 || MI(mm1) == 4
        x1 = (MY(mm1) - b1)/m;
        x2 = (MY(mm1) - b2)/m;
        if MX(mm1) >= x2 && MX(mm1) <= x1
          %  MI(mm1) = 11;
        end
    end       
    
end


%--------------------------------------------------------
% Initial temperature structure
%--------------------------------------------------------

mm1 = MY >= cont_lev;
MTK(~mm1) = ttop;
MTK(mm1) = interp1(SE.z,SE.T,MY(mm1)-cont_lev);
    
% -------------------------------------------------------------------------
%% Initial elevation for topography profile
% -------------------------------------------------------------------------

% Intial elevation at bottom of sticky layer
gridt(2,:) = sticky_layer;

% Set initial continent elevation
k = gridt(1,:) >= lx1 & gridt(1,:) <= lx2;
gridt(2,k) = sticky_layer - lh;

% Save initial topography
% topotime(1,1) = 0;
% topohigh(1,:) = gridt(2,:);
% if ~isempty(water_lev)
%     topo_lev=zeros(ntimestep+1,tnum);
%     topowater(1,1) = water_lev;
% end