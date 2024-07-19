function Prasanna_RunThermoMec2D(gridmodel,inputmodel,modelbc,outpath)

% =======================================================================
% RunThermoMec2D
% =======================================================================
%
% Version 2.1 Feb 7, 2017
% Rob Moucha, Syracuse University
% Siobhan Campbell, Syracuse University
% Prasanna Gunawardana, Syracuse University
%
% Added seperate input script file for boundary conditions
% - Bug fixes, melt and fluid weakening corrected
%
% Version 2.0 Oct. 2016
% Robert Moucha, Syracuse Unviersity
% Siobhan Campbell, Syracuse University
% 
% Fixed bug with boundary marker integration
% New Boundary Condition Model File
% 
% Thermomechanical visco-elasto-plastic numerical model with
% erosion/sedimentation processes based on Gerya (2010).
% 
%
% Version 1.0, Sep. 2014
% Robert Moucha, Syracuse University
% Peter Nelson,  University of Texas at Austin
%
% External functions used:
%   <Model Setup File>
%   ThermoMec2D_Fast.m
%   Fast_Stokes_Continuity_solver_sandbox.m,
%   Fast_Temperature_solver_grid.m
%   marker_rheology.m
%   markers2grid.m
%   movemarkers.m
%   
%   
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Solve Primitive Variable Stokes flow on a staggered grid with
%
%     vx       vx       vx
%
% vy  T---vy---T---vy---T   vy20
%     |        |        |
%     vx   P   vx   P   vx
%     |        |        |
% vy  T---vy---T---vy---T   vy
%     |        |        |
%     vx   P   vx   P   vx
%     |        |        |
% vy  T---vy---T---vy---T   vy
%
%     vx       vx       vx
%
% Lines show basic grid, ghost nodes shown outside the basic grid are used
% for boundary conditions.
% -------------------------------------------------------------------------

% Amount of timesteps
stepmax = 2100; %% Exceso, se define con el tiempo (model time)1300

% Desired model time (yrs);
modeltime = 51.0e6;

% Frequency of output / Big and small segment
nfreq_big = 20;
nfreq_small = 100000;

% limits of small output box:
markx=[1200 2000]*1000;
marky=[0 70]*1000;

if ~exist(outpath,'dir')
    mkdir(outpath)
end

% Parameter files 
% Save copies of the input model and grid model in the output path
infile_copy = [outpath,'/input_model.m'];
copyfile([inputmodel,'.m'],infile_copy);
grid_copy = [outpath,'/grid_model.m'];
copyfile([gridmodel,'.m'],grid_copy);

% Get the full path to the output directory
curdir = pwd;
cd(outpath)
outpath = pwd;
cd(curdir);

% Remove horizontal boundary conditions after period of time remove_time%
remove_time =3.5;   % Ma, set to false to turn-off% Strechtching try 3.5Ma%false
%remove_time1 = 3.0;
%remove_time2 = 5.0;
%%%%%%%%%%%%%%%%%%%%%%%%%
% Add plume head at a point in time, give plume parameters in input_model script
%add_plume = true;  % Ma, set to false to turn-off%Kepp this as true and change the add_plume in the input model

% Maximal temperature change, allowed for one timestep (K)
% NOTE: If tempmax = 0 ==> Stokes Flow only
tempmax = 35;        % Maximum temperature change per timestep, if exceeded heat
                     % heat equation time step will be reduced
abstempmax = 2000;   % Unphysical temperature change, something very wrong

if tempmax==0
    warning('Stokes Flow Only Simulation')
end

% Shear heating on(1)/off(0)
frictyn = 1;

% Adiabatic heating on(1)/off(0)
adiabyn = 1;

% Restart Flag, must be false otherwise it will not work
restart = false;

% Melting flag, melt on(true)/off(false)
melting = true;

% Pore fluid and melt pressure factors (Sizova et al., 2010; Gerya and Meilick, 2011)
% The pore fluid pressure P_fluid reduces the yield strength of fluid
% bearing rock. For dry rocks lambda_fluid = 1. Hydrostatic gradiant in the
% upper crust is lamda_fluid = 0.6. Likewise, lamda_melt is due to weaking
% of rocks due to ascending melts.
lambda_fluid = 0.001; % 1 - P_fluid/P_solid
lambda_melt = 0.001;  % 1 - P_melt/P_solid

% Radius of ascenet channel/dike (m) 
dike_dx = 50;

% Partially Melt Fractions
% Melt extraction threshold fraction
meltmax = 0.04;
% Non-extractable amount of melt
meltmin = 0.02;
% Partially molten rock viscosity
etamelt=1e+17; % Pa s

% -------------------------------------------------------------------------
%% Constants and conversion factors
% -------------------------------------------------------------------------

% Acceleration of Gravity (m/s^2)
gx=0;
gy=9.81;

% Gas constant (J/mol/K)
RGAS=8.314;

% Conversion year to second
yr2sec = 365.25*24*3600;

% Absolute zero (deg-C)
tabsolute = -273;

% Adiabatic temperature gradient K/km
tgrad = 0.5;
tgrad = tgrad/1000; % K/m

% -------------------------------------------------------------------------
%% Initial Computational Grid and Marker Settings
% -------------------------------------------------------------------------
% Run input grid file to set grid parameters
run(gridmodel)

% Total number of material markers per cell
mxcell = 5;
mycell = 5;

mxnum = mxcell*(xnum-1);  % total number of desired markers in horizontal direction
mynum = mycell*(ynum-1);  % total number of desired markers in vertical direction
xsize = xsize0;           %% tamaño físico del modelo   en metros, tomado de gridmodel
ysize = ysize0;

% -------------------------------------------------------------------------
% Initial temperature at the top, and bottom of the model based on2
% prescribed mantle potential temperature (oC)
% -------------------------------------------------------------------------
ttop = 0;          % (deg-C)
tpotential = 1343; % (deg-C)%1343

% -------------------------------------------------------------------------
%% Viscosity and stress limits for rocks
% -------------------------------------------------------------------------
viscoelastic = true;     % True or False to use viscoelastic formulation  KEEP TRUE (false still needs to be fully implemented)
etamin = 1.0e+18;        % Lower limit (Pa)
etamax = 1.0e+25;        % Upper limit (Pa)
stressmin = 1.0e+4;      % Lower stress limit for power law (Pa)

timemax = 1.0e+4;%3.5e+4;        % see Viscoelastic timestep (5,000 yr) p. 183 Gerya

timemax = timemax*yr2sec; % Viscoelastic timestep (s) p. 183 Gerya
% !!! If you see an ossicilation of max velocity, reduce time step see p. 275 !!!

% -------------------------------------------------------------------------
%% Marker Motion, Subgrid Diffusion
% -------------------------------------------------------------------------
% Maximal marker displacement step, number of gridsteps
markmax=1;%0.5
% Moving Markers:
% 0 = not moving at all
% 1 = simple advection
% 4 = 4-th order in space  Runge-Kutta
markmove=4;
% Velocity calculation
% 0 = by Solving momentum and continuity equations
% 1 = solid body rotation
% Velocity calculation
% 0 = by Solving momentum and continuity equations
% 1 = solid body rotation
movemod=0;
% Numerical Subgrid stress diffusion coefficient
dsubgrids=1;
% Numerical Subgrid temperature diffusion coefficient
dsubgridt=1;

% Recycle markers across boundary
markerreflect = true;

% Fill empty cells with markers of the type
fillemptycells = false;
filllastrow = false;
newmarkertype = 10; % Asthenospheric mantle
addmarkers = false;

% ------------------------------------------------------------------------
%% Topography and Surface Processes
% -------------------------------------------------------------------------

% Topography model size in horizontal direction
tsize = xsize0;
% Defining topography model resolution
tstp = 500; % (m)
tnum = ceil(xsize0/tstp + 1);

% Grid for topography profile: 1 => x, 2 => y , 4 => vx, 5 => vy, 
%                              3,6 = >contianers for surface processes
gridt =zeros(6,tnum);
gridt(1,:) = linspace(0,xsize0,tnum);
tstp = gridt(1,2) - gridt(1,1); % True resolution

% Define erosion rate (dYt/dt) 
dYtdt = 0.5; % (mm/yr)
dYtdt = dYtdt/1000/yr2sec; % (m/s)

% Define max elevation (d2Yt) 
d2Yt = 10000; % (m)

% Define transport lengthscale (dx)
dx = 100000; % (m)

% Topography diffusion koefficient Ks (m^2/s)
% dYt/dt=Ks*d2Yt/dx^2
Ks = dYtdt*dx^2/d2Yt;

% Allocate topography containers
% topotime = zeros(1000,1);
% topohigh = zeros(1000,tnum);
% topowater = zeros(1000,1);

% -------------------------------------------------------------------------
%% Construct Computational Grid & Populate with blank Markers
% -------------------------------------------------------------------------

% Assamble non uniform grid parameters
G.xnum = xnum; G.ynum = ynum;
G.xsize = xsize; G.ysize = ysize;
G.bx = bx; G.by = by; % Resolution of fine region (m)
G.wx = wx; G.wy = wy; G.f = f;  % width of fine uniform grid (m)

% Construct the grid
G = construct_grid(G);
grid_type = G.grid_type;
gridx = G.x'; gridy = G.y';
gridx1 = gridx(1);

% Defining intial position of markers
M.nx = mxnum; M.ny = mynum;
M.mxcell = mxcell;
M.mycell = mycell;
M = generate_markers(M,G);
MX = M.x; MY = M.y;
clear M G

G.grid_type = grid_type;
G.xsize = gridx(end);
G.ysize = gridy(end);
G.xnum = xnum;
G.ynum = ynum;
G.bx = bx; G.by = by; G.wx = wx; G.wy = wy; G.f = f;
G.gridx = gridx;
G.gridy = gridy;

% Creating markers array
marknum = mxnum*mynum;

MTK=zeros(marknum,1);   % Temperature, K
MI=zeros(marknum,1);    % Type
MXN=zeros(marknum,1);   % Horizontal index
MYN=zeros(marknum,1);   % Vertical index
MCXN=zeros(marknum,1);  % Horizontal central index
MCYN=zeros(marknum,1);  % Vertical central index
MSXX=zeros(marknum,1);  % SIGMAxx - deviatoric normal stress, Pa
MSXY=zeros(marknum,1);  % SIGMAyy - shear stress, Pa
META=zeros(marknum,1);  % viscosity, Pa s
MEXX=zeros(marknum,1);  % EPSILONxx - normal strain rate, 1/s
MEXY=zeros(marknum,1);  % EPSILONyy - shear strain rate, 1/s
MPR=zeros(marknum,1);   % Pressure, Pa
MGII=zeros(marknum,1);  % Accumulated strain
MRAT=ones(marknum,1);   % EiiMarker/EiiGrid Ratio
MXM=zeros(marknum,1);   % Cummulative Melt Fraction
MEXTC=zeros(marknum,1); % Cummulative Extracted Melt Fraction
MEXT=zeros(marknum,1);  % Extracted Melt Fraction
MLAMBDA=ones(marknum,1); % Weakening effect of melt/fluid

% -------------------------------------------------------------------------
%% Define Material Properties & Build Geometry
% -------------------------------------------------------------------------

% Run input parameter file
run(inputmodel)

% Run to generate boundary conditions
run(modelbc)

%% Run ThermoMec2D
ThermoMec2D_Fast
