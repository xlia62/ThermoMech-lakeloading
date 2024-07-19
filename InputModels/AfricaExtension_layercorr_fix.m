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

AfricaMaterials_layercorr
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

% Turn off melting
melting = false;

% Sticky Layer thickness water or air (m)
sticky_layer = 20000;

% Water Level Tracking (m) set to empty if no water
% Water depth = sticky_layer - water_lev
water_lev = [];%10750;%20200;

% water_lev will be adjusted according to max water depth (m)
% Water depth will be set to 
max_water_depth = 1000;

% Sediment Thickness (m)
sed_thick = 0;

% Continent level (m), use this to set continent above or at sea level
cont_lev = sticky_layer;

% Upper Continental Crust thickness (m)
upper_cont_thick = 11000 + sed_thick;  % Thickness includes sediment layer

% How many layers in upper crust to alternate for visualizing deformation
num_crust = upper_cont_thick/1000;

% Lower Continental Crust thickness (m)
lower_cont_thick = 25000;

% Crustal thickness (with sediments) (m)
crust_thick = upper_cont_thick + lower_cont_thick;

% Mantle Lithosphere (m) 
mantle_lith = 105000;

% Lithosphere (LAB) (m) % Specifically, the Thermal lithosphere depth
LAB_depth = mantle_lith + crust_thick;

% Weak Zones

% Weak seed
xweak = 1000;
yweak = 1000;
x0 = xsize/2;
y0 = 5000 + sticky_layer;
x1 = x0 - xweak/2;
x2 = x0 + xweak/2;
y1 = y0 - yweak/2;
y2 = y0 + yweak/2;

% Size and position of half-disk weak zone (m)
%rweak = 5000;
%xweak = xsize/2;
%yweak = upper_cont_thick  + sticky_layer; % brittle region of mantle lithosphere (model specific)

% Fault
% dip_angle = 60;    %(deg) <-- change to 0 if want single spot of weakness
% fault_width = 1000; % (m) <-- for fault = 500; spot == 1000
x0 = xsize/2; 
y0 = sticky_layer; % Change to center of UC if want single spot of weakness
% y0 = sticky_layer + upper_cont_thick/2 -500;
% m = -tand(dip_angle);
% b1 = y0 - m*x0;  % Left bound intercept
% b2 = b1 + m*fault_width;  % Right bound intercept

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
global dtpotential
%dtpotential = 0;%100;

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

global tmoho
switch n
    
    case 0
        
        % Moho temperature 700oc
        S0 = tmoho;
        %S0 = 600;          % Moho Temperature (oC)
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
        Zrad = 10000; % Constant Radioactive zone % Use for Constant heat production rate
        qm = 0.03; % Heat flow of the continent - ??
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
        MI(mm1)= 3;
    end
    
    % Continental crust
    % Bottom of upper crust
    bot_crust1 = bot_sed + upper_cont_thick;
    if(MY(mm1)>bot_sed && MY(mm1)<=(bot_crust1-upper_cont_thick/2-500)...
            || MY(mm1)>(bot_crust1-upper_cont_thick/2+500) && MY(mm1)<=bot_crust1) 
        
        MI(mm1) = 6;
        
        % Layered structure in the sand
        m2 = double(int16(MY(mm1)/(upper_cont_thick/num_crust)-0.5)) + 1;
        m2 = m2 - double(int16(m2/2-0.5))*2;
        if(m2==0)
            MI(mm1)= 4;
        end
    end
    
    % Tracking layer in Upper Crust
    if (MY(mm1)>=bot_crust1-upper_cont_thick/2-500 && MY(mm1)<=bot_crust1-upper_cont_thick/2+500)
        
        MI(mm1) = 5;
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
    
%     % Fault zone
%     if MI(mm1) == 6 || MI(mm1) == 4 || MI(mm1) == 5
%         x1 = (MY(mm1) - b1)/m;
%         x2 = (MY(mm1) - b2)/m;
%         if MX(mm1) >= x2 && MX(mm1) <= x1
%             MI(mm1) = 11;
%         end
%     end     
    
    % Weak spot <-- use if seeding a spot rather than fault
    if MY(mm1) >= y1 && MY(mm1) <= y2
        if MX(mm1) >= x1 && MX(mm1) <= x2
            MI(mm1) = 11;
        end
    end
end


%--------------------------------------------------------
% Initial temperature structure
%--------------------------------------------------------

mm1 = MY >= cont_lev;
MTK(~mm1) = ttop;
MTK(mm1) = interp1(SE.z,SE.T,MY(mm1)-cont_lev);
    
% ------------------------------------------------------------------------
%% Topography and Surface Processes
% -------------------------------------------------------------------------

% Use a surface process timestep <5000 yrs.
LEMpar.dt_max = 2000; %yrs

% Apply surface processes model true/false
LEMpar.apply_surfacepro = true;

% Define erosional parameters
LEMpar.Kf = 2e-5; % bedrock river incision rate meters^(1-2m)/yr
LEMpar.Kd = 0.1;  % bedrock transport coefficient (diffusivity) meters^2/yr
LEMpar.Fd = 0;  % Effeciency of fluvial sediment deposition. 
                  % You can turn off fluvial deposition by setting to zero.        
LEMpar.Ld = false; % Deposit sediment in lakes true/false. 

% Desired topography model resolution
tstp = 250; % (m)

% Lake depth is set-up as a global paremeter
global water_depth
%water_depth = 1200;

% Topography model size in horizontal direction, adjust resolution to fit
tsize = xsize0;
tnum = ceil(xsize0/tstp + 1);

% Grid for topography profile: 1 => x, 2 => y , 4 => vx, 5 => vy,
%                              3,6 = >contianers for surface processes
%                              7 no surface process topography
gridt =zeros(7,tnum);
gridt(1,:) = linspace(0,xsize0,tnum);

% True topography resolution
tstp = gridt(1,2) - gridt(1,1); 

% Intial elevation at bottom of sticky layer
gridt(2,:) = sticky_layer;

% Set-up eroisonal surrace and stratigraphy horizons

% Erosional surface, tracks surface as if no erosion/deposition
E = gridt(2,:);

% Age interval
dt_horizon = 1e6;  % yr
dt_horizon = dt_horizon*yr2sec;  % sec.

% Number of horizons
n_horizon = floor(Par.modeltime/dt_horizon)+1;  

% The horizon array, 1st row corresponds to t = dt_horizon.
% Horizon is initialized in input model
i_horizon = 1;
H = zeros(n_horizon,tnum);

% The first horizon is the basement.
H(1,:) = gridt(2,:);
t_horizon = dt_horizon;  %Set next horizon time.
