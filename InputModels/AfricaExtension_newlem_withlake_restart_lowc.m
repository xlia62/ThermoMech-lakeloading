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
TheunissenMaterials_withlake_lowC

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
step2restart = 110; 
load(['../output/AfricaModels2022/',...
'Lake_drop_restart_thermodynamic_reference_lowc/markers_',num2str(step2restart),'.mat'])
load(['../output/AfricaModels2022/',...
'Lake_drop_restart_thermodynamic_reference_lowc/grids_',num2str(step2restart),'.mat'])
load('../output/AfricaModels2022/Lake_drop_restart_thermodynamic_reference_lowc/bounds.mat')

% load('/home/lxue07/Documents/ThermoMech2d/ForLiang_ThermoMech2D/ForLiang_ThermoMech2D/output/AfricaModels2022/Lake_drop_restart_test2/markers_',num2str(step2restart),'.mat')
% load('/home/lxue07/Documents/ThermoMech2d/ForLiang_ThermoMech2D/ForLiang_ThermoMech2D/output/AfricaModels2022/Lake_drop_restart_test2/grids_30.mat')
% load('/home/lxue07/Documents/ThermoMech2d/ForLiang_ThermoMech2D/ForLiang_ThermoMech2D/output/AfricaModels2022/Lake_drop_restart_test2/bounds.mat')
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
water_lev = 20000; %[];%10750;%20200;
water_lev_bot = 20800; % min of lake level, 800 m water level change, was 24000
period = 400; %kyr 
numer_of_p = 2; % number of period in the model 
buffer_time = 0.2; %Myr time after restart, and before lake level changes

%% Topography and Surface Processes
% True topography resolution
% Use a surface process timestep <5000 yrs.
LEMpar.dt_max = 2000; %yrs

% Apply surface processes model true/false
LEMpar.apply_surfacepro = false;

% Define erosional parameters
LEMpar.Kf = 5e-6; % bedrock river incision rate meters^(1-2m)/yr
LEMpar.Kd = 1;    % bedrock transport coefficient (diffusivity) meters^2/yr
LEMpar.Ks = 500;  % marine sediment transport coefficient
LEMpar.Fd = 1;    % Effeciency of fluvial sediment deposition. 
                  % You can turn off fluvial deposition by setting to zero.  
LEMpar.sea_level = water_lev;  % Define sea_level for surface processes and marine deposition                
LEMpar.Ld = true; % Deposit sediment in lakes true/false. 

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

%% Initial Temperature parameters
%--------------------------------------------------------------------------

% Plume Excess Temperature [K]
add_plume = false;  % Ma,  set to false







