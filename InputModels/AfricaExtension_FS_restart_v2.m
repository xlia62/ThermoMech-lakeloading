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

%AfricaMaterials_layercorr
TheunissenMaterials_withlake

% The materials have 2 different markers for upper continental crust to
% visualize deformation

% Define plastic healing characterstic time, plastic_strain = plastic_strain/(1+timestep/tau_plastic);
tau_plastic = 15e6; % years (10^10-10^12), or set to very high for no plastic healing
tau_plastic = tau_plastic*yr2sec;

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

RxType =  'thermodynamic'; %; was false
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
%% load marker and grid
% -------------------------------------------------------------------------
step2restart = 69; 
load(['../output/AfricaModels2022/',...
'Lake_drop_model_reference_fine/markers_',num2str(step2restart),'.mat'])
load(['../output/AfricaModels2022/',...
'Lake_drop_model_reference_fine/grids_',num2str(step2restart),'.mat'])
load('../output/AfricaModels2022/Lake_drop_model_reference_fine/bounds.mat')

startstep = ntimestep;

% time of lake level starting to change 
water_lev_beg_t = timesum; % second 

%% Defining lithological model constants
% -------------------------------------------------------------------------

% Sticky Layer thickness water or air (m)
sticky_layer = 20000;

% Sediment Thickness (m)
sed_thick = 0;

% Continent level (m), use this to set continent above or at sea level
cont_lev = sticky_layer;

% Upper Continental Crust thickness (m)
upper_cont_thick = 15000 + sed_thick;  % Thickness includes sediment layer

% Lower Continental Crust thickness (m)
lower_cont_thick = 25000;

% Crustal thickness (with sediments) (m)
crust_thick = upper_cont_thick + lower_cont_thick;

% How many layers in upper crust to alternate for visualizing deformation
num_crust = fix(crust_thick/500);

% Mantle Lithosphere (m) 
mantle_lith = 60000;

% Lithosphere (LAB) (m) % Specifically, the Thermal lithosphere depth
LAB_depth = mantle_lith + crust_thick;

%% Generate water level change @ Liang
%--------------------------------------------------------------------------
% water level change based on Milankovitch cycle
% eccentricity period of 400 kyr
global lake_level_change
% lake_level_change = true;

% Sea level Tracking (m) set to empty if no water
% Water depth = sticky_layer - water_lev
water_lev = 22000; %[];%10750;%20200; max of lake level 
water_lev_bot = 22600; % min of lake level, 800 m water level change, was 24000
ini_water_lev =  water_lev;
period = 100; %kyr 
numer_of_p = 2; % number of period in the model 
buffer_time = 0.075; %Myr time after restart, and before lake level changes

%% Topography and Surface Processes
% True topography resolution
% Use a surface process timestep <5000 yrs.
LEMpar.dt_max = 2000; %yrs

% Offset used to map from ThermMech grid (positive down) to sea-level at 0 and positive up.
LEMpar.offset = sticky_layer;
LEMpar.sea_level = LEMpar.offset-water_lev;  % Define sea_level for surface processes and marine deposition, if water_lev is empty, it will be empty 
LEMPar.max_sea_level = LEMpar.sea_level;     % Define the absolut max of sea level, important for dynamic changes driven by max_depth.
LEMpar.max_depth = 500;  % Define maximum depth (m), this will cause the LEMpar.sea_level and water_lev to change accordingly 
                         % set to very high if no control is desired                        
                         
% Apply surface processes model true/false
LEMpar.apply_surfacepro = 'none';

% Define erosional parameters
LEMpar.Kf = 0;  %5e-6; % bedrock river incision rate meters^(1-2m)/yr
LEMpar.Kd = 0; %1;    % bedrock transport coefficient (diffusivity) meters^2/yr
LEMpar.Ks = 0; %20;  % marine sediment transport coefficient
LEMpar.Fd = 0; %2;    % Effeciency of fluvial sediment deposition. 
                  % You can turn off fluvial deposition by setting to zero.  
% LEMpar.sea_level = water_lev;  % Define sea_level for surface processes and marine deposition

LEMpar.Ld = false; % Deposit sediment in lakes true/false. 
LEMpar.Lcrit = 2000; % Marine Deposition critical distance

% Additional FastScape erosional parameters
LEMpar.m = 0.5;
LEMpar.n = 1.0;
LEMpar.kfsed = -1;
LEMpar.kdsed = -1;
LEMpar.expp = -2;
LEMpar.p = 1;

% Desired topography model resolution
tstp = 250; % (m)

% Lake depth is set-up as a global paremeter
global water_depth
if isempty(water_depth)
    water_depth = 0;
end


%  Defining new gridline positions for irregular basic grid
xsize0 = G.xsize;
ysize0 = G.ysize;

% Topography model size in horizontal direction, adjust resolution to fit
% tsize = xsize0;
% tnum = ceil(xsize0/tstp + 1);
tsize = gridt(1, 1) - gridt(1, end);
tnum = length(gridt);

% Age interval
dt_horizon = 1e6;  % yr
dt_horizon = dt_horizon*yr2sec;  % sec.
t_horizon = dt_horizon;  %Set next horizon time.

% Horizon is initialized in input model
i_horizon = 1;


%% 
% Use a surface process timestep <5000 yrs.
% % % LEMpar.dt_max = 2000; %yrs
% % % 
% % % % Offset used to map from ThermMech grid (positive down) to sea-level at 0 and positive up.
% % % LEMpar.offset = sticky_layer;
% % % LEMpar.sea_level = LEMpar.offset-water_lev;  % Define sea_level for surface processes and marine deposition, if water_lev is empty, it will be empty 
% % % LEMPar.max_sea_level = LEMpar.sea_level;     % Define the absolut max of sea level, important for dynamic changes driven by max_depth.
% % % LEMpar.max_depth = 500;  % Define maximum depth (m), this will cause the LEMpar.sea_level and water_lev to change accordingly 
% % %                          % set to very high if no control is desired                        
% % % 
% % % % Apply surface processes model: 'none','line','fastscape'
% % % % Make sure you edit the fastscape setup at the bottom of this file if
% % % % using 'fastscape'
% % % LEMpar.apply_surfacepro =  'none';
% % % 
% % % % Define erosional parameters
% % % LEMpar.Kf = 5e-6; % bedrock river incision rate meters^(1-2m)/yr
% % % LEMpar.Kd = 1;    % bedrock transport coefficient (diffusivity) meters^2/yr
% % % LEMpar.Ks = 20;   % marine sediment transport coefficient
% % % LEMpar.Fd = 2;    % Effeciency of fluvial sediment deposition (G in fastscape)
% % %                   % You can turn off fluvial deposition by setting to zero.  
% % %              
% % % LEMpar.Ld = false ; % Deposit sediment in lakes true/false. 
% % % LEMpar.Lcrit = 2000; % Marine Deposition critical distance
% % % 
% % % % Additional FastScape erosional parameters
% % % LEMpar.m = 0.5;
% % % LEMpar.n = 1.0;
% % % LEMpar.kfsed = -1;
% % % LEMpar.kdsed = -1;
% % % LEMpar.expp = -2;
% % % LEMpar.p = 1;
% % % 
% % % % Desired topography model resolution
% % % tstp = 250; % (m)



% % Topography model size in horizontal direction, adjust resolution to fit
% tsize = xsize0;
% tnum = ceil(xsize0/tstp + 1);
% 
% % Grid for topography profile: 1 => x, 2 => y , 4 => vx, 5 => vy,
% %                              3,6 = >contianers for surface processes
% %                              7 no surface process topography
% gridt =zeros(7,tnum);
% gridt(1,:) = linspace(0,xsize0,tnum);
% 
% % True topography resolution
% tstp = gridt(1,2) - gridt(1,1); 



%% Initial Temperature parameters
%--------------------------------------------------------------------------

% Plume Excess Temperature [K]
add_plume = false;  % Ma,  set to false







