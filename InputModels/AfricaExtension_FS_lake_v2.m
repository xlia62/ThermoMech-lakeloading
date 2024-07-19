% Define Material Properties, Build Geometry, Construct Boundary Conditions
%
% Edit this file to change:
%
%         material properties
%         model geometry
%         boundary conditions
%         landscape evolution parameters

% -------------------------------------------------------------------------
%% Material properties
% -------------------------------------------------------------------------

TheunissenMaterials_withlake

% Define plastic healing characterstic time, plastic_strain = plastic_strain/(1+timestep/tau_plastic);
tau_plastic = 15e6; % years (10^10-10^12), or set to very high for no plastic healing
tau_plastic = tau_plastic*yr2sec;

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

RxType = 'thermodynamic'; %'thermodynamic';
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

% Sea level Tracking (m) set to empty if no water
% sea_level = sticky_layer - water_lev
water_lev = 20000;%10750;%20200;

% Sediment Thickness (m)
sed_thick = 0;

% Continent level (m), use this to set continent above or at sea level
cont_lev = sticky_layer;

% Upper Continental Crust thickness (m)
upper_cont_thick = 20000 + sed_thick;  % Thickness includes sediment layer

% Lower Continental Crust thickness (m)
lower_cont_thick = 20000;

% Crustal thickness (with sediments) (m)
crust_thick = upper_cont_thick + lower_cont_thick;

% How many layers in upper crust to alternate for visualizing deformation
num_crust = fix(crust_thick/500);

% Mantle Lithosphere (m) 
global mantle_lith
mantle_lith = 60000;
if  isempty (mantle_lith)
    mantle_lith = 60e3;
end

% Lithosphere (LAB) (m) % Specifically, the Thermal lithosphere depth
LAB_depth = mantle_lith + crust_thick;

%--------------------------------------------------------------------------
%% Generate water level change @ Liang
%--------------------------------------------------------------------------
global lake_level_change 
% true if water level change, here we use density change 
% to estimate the water level changes.

%lake_level_change = true;
water_density_ini = MRHO(2,1);

% ice box model variables
% Ice height/ rift lake depth(m)
% ice_height = 2000;

% water wedge model variables
% % the x limits of ice
wedge_center= 130e3;
xl_1 = 110e3;
xl_2 = 150e3;
% wedge_center= 175e3;
% xl_1 = 150e3;
% xl_2 = 200e3;
% basin depth at center, m
basin_d = 2.5e3;  


% ice_decrease_t_ini = 23; % Ma time when the lake level starts to change
% ice_decrease_t_end = 23.5; % Ma time when the lake level stops to change
% 
% % the density change would be 
% % MRHO(2,1) = ice_density_ini * (1 + t_since_change * ice_level_change_factor);
% ice_level_change_factor = -2e-6; % percentage/ yr +: increase, -: decrease


%--------------------------------------------------------------------------
%% Setup initial Weak Zones
%--------------------------------------------------------------------------
% Size and position of half-disk weak zone (m)
% rweak = 5000;
%  xweak = xsize/2;
%  yweak = crust_thick  + sticky_layer - rweak; % brittle region of mantle lithosphere (model specific)

% Seed type: Weak seed = 1, fault dip = 2, tapered random weakness = 3 
seed_type = 2;
thermal_weakness = false; %add a thermal perturbation at base of lithosphere

if seed_type == 1
    xweak = 2000;
    yweak = 2000;
    x0 = G.xsize/2;
    y0 = 6000 + sticky_layer - yweak/2;
    x1 = x0 - xweak/2;
    x2 = x0 + xweak/2;
    y1 = y0 - yweak/2;
    y2 = y0 + yweak/2;
elseif seed_type == 2
    % Fault
    dip_angle = 60;    %(deg)
    fault_width = 1000; %(m)
    x0 = G.xsize/2;
    y0 = sticky_layer;
    m = -tand(dip_angle);
    fault_bot = crust_thick + sticky_layer + 5e3; % depth of fault penetrate crust 
    b1 = y0 - m*x0;  % Left bound intercept
    b2 = b1 + m*fault_width;  % Right bound intercept
elseif seed_type == 3 
    % width of weakness zone in m
    x0 = G.xsize/2;    %center of weak zone
    zweak = 100000;  %width of weak zone
    x1 = x0 - zweak/2;
    x2 = x0 + zweak/2;
    x_taper = linspace(x1,x2,201);
    y_taper = cosine_tapered_win(length(x_taper),0.25);
end

if thermal_weakness
    rweak = 20000;  %radius of weakzone (m)
    % Weak zone temperature perturbation [k];
    Tweak = 20;
end

%--------------------------------------------------------------------------
%% Initial Temperature parameters
%--------------------------------------------------------------------------

% Plume Excess Temperature [K]
add_plume = false;  % Ma,  set to false

% Thermal diffusivity of the mantle, m^2/s
kappa = 1e-6;    

% Increase default mantle potential temperature oC
global dtpotential
if isempty(dtpotential)
    dtpotential = 0;
end

% Convert deg-C to K
ttop = ttop - tabsolute;
tpotential = tpotential + dtpotential - tabsolute;  % defined at RunThermoMec2D

% Temperature at bottom of model space
tbottom = 1200 - tabsolute;  % Will be updated below after the

%--------------------------------------------------------------------------
%% Generate initial geotherm and continental strength envelopes
%--------------------------------------------------------------------------

% Continental rifting only -- needs to be modified for general case
% 
% Assemble The geotherm (Temperature Profile to LAB)
zmax = G.ysize-sticky_layer;    % Depth of geotherm for strength envelop and lithosphere
SE.z = linspace(0,zmax,500)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assemble the geotherm
% n= 0 for piece-wise linear, 1 for exponential heat production 2 for constant heat production
n = 0;

Zl = LAB_depth;    % Mechanical Lithosphere Depth (km)
Ti = ttop + tabsolute; % Ti - Surface Temperature (oC)
Tl = 1200;             % Tl = Temperature of the base of mechanical Mantle Lithosphere (oC)

% Moho Temperature (oC)
global tmoho
if isempty(tmoho)
    tmoho = 600;
end
switch n
    
    case 0
        
        % Moho temperature
        S0 = tmoho;
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
SE.eii = 1e-15;
% Continental Lithosphere
SE = strengthprofile(sed_thick,MRHO(3,:),MFLOW(3,:),MPL(3,:),...
    upper_cont_thick-sed_thick,MRHO(6,:),MFLOW(6,:),MPL(6,:),...
    lower_cont_thick,MRHO(6,:),MFLOW(6,:),MPL(6,:),...  
    mantle_lith,MRHO(9,:),MFLOW(9,:),MPL(9,:),...
    zmax,MRHO(10,:),MFLOW(10,:),MPL(10,:),SE);
%lower_cont_thick,MRHO(7,:),MFLOW(7,:),MPL(7,:),
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
    % !!! lower and upper crust are the same !!!
    
    crust_thick = upper_cont_thick + lower_cont_thick;
    bot_crust1 = bot_sed + upper_cont_thick;
    bot_crust2 = bot_crust1 + lower_cont_thick;
    
    if(MY(mm1)>bot_sed && MY(mm1)<=bot_crust2)
        
        MI(mm1) = 6;
        % Layered structure in the sand
        m2 = double(int16(MY(mm1)/(crust_thick/num_crust)-0.5)) + 1;
        m2 = m2 - double(int16(m2/2-0.5))*2;
        if(m2==0)
            MI(mm1)= 4; % @Liang 
        end
    end
    
    % Strong Continetal Mantle Lithosphere
    if(MY(mm1)>bot_crust2 && MY(mm1)<= LAB_depth)
        MI(mm1) = 9;
    end
   
    % Lower mantle
    if(MY(mm1) > 660e3+sticky_layer)
        MI(mm1) = 15;
    end
    
    %--------------------------------------------------------
    %  Weak seed
    %--------------------------------------------------------
    
    % Thermal perturbation in the base of the lithosphere
    % (to start extension in the middle of the model)
    
    
    %     if MY(mm1) >= y1 && MY(mm1) <= y2
    %         if MX(mm1) >= x1 && MX(mm1) <= x2
    %             MI(mm1) = 11;
    %         end
    %     end
    
    % Weak spot <-- use if seeding a spot rather than fault
    if seed_type == 1
        if MY(mm1) >= y1 && MY(mm1) <= y2
            if MX(mm1) >= x1 && MX(mm1) <= x2
                MI(mm1) = 11;
            end
        end
    elseif seed_type == 2
        % Fault zone
        if MY(mm1)<=fault_bot && (MI(mm1) == 6 || MI(mm1) == 4|| MI(mm1) == 9|| MI(mm1) == 10)
            x1 = (MY(mm1) - b1)/m;
            x2 = (MY(mm1) - b2)/m;
            if MX(mm1) >= x2 && MX(mm1) <= x1
                MI(mm1) = 11;
            end
        end
    elseif seed_type == 3
        if MY(mm1) > sticky_layer && MY(mm1)<=LAB_depth
            if MX(mm1) >= x1 && MX(mm1) <= x2
               mgii = 0.15 + randn/20;  
               if mgii>0
                   if mgii > 0.3
                       mgii = 0.3;
                   end
                   MGII(mm1) = mgii*interp1(x_taper,y_taper,MX(mm1));
               end
            end
        end
    end
    
%-------------------------------------------------------------
%%  @Liang Ice/Rift lakes, define after sticky are defined
%-------------------------------------------------------------
    % wedge like basin 
    % tan dipping 
    tand_d = basin_d/((xl_2-xl_1)/2);  
    
    % left side of the wedge    
    if (MY(mm1)>=sticky_layer && MY(mm1) <= sticky_layer + basin_d && MX(mm1)>= xl_1 && MX(mm1)<= wedge_center)
        y_d = (MX(mm1)-xl_1)*tand_d + sticky_layer;
        if (MY(mm1) <= y_d)
            MI(mm1) = 2;
        end
    end    
    
    % right side of the wedge
    if (MY(mm1)>=sticky_layer && MY(mm1) <= sticky_layer + basin_d && MX(mm1)> wedge_center && MX(mm1)<= xl_2)
        y_d = (xl_2 - MX(mm1))*tand_d + sticky_layer;
        if (MY(mm1) <= y_d)
            MI(mm1) = 2;
        end
    end   
    
    
end

%--------------------------------------------------------
% Initial temperature structure
%--------------------------------------------------------

mm1 = MY >= cont_lev;
MTK(~mm1) = ttop;
MTK(mm1) = interp1(SE.z,SE.T,MY(mm1)-cont_lev);

if thermal_weakness
    for mm1 = 1:marknum
        % Thermal perturbation in the base of the lithosphere
        % (to start extension in the middle of the model)
        dx = MX(mm1) - G.xsize/2;
        dy = MY(mm1) - LAB_depth;
        dr = (dx^2+dy^2)^0.5;
        if (MY(mm1) <= LAB_depth && dr <= rweak)
            MTK(mm1) = MTK(mm1) + Tweak*(1-dr/rweak);
        end
    end
end
%

% ------------------------------------------------------------------------
%% Topography and Surface Processes Parameters
% -------------------------------------------------------------------------

% Use a surface process timestep <5000 yrs.
LEMpar.dt_max = 2000; %yrs

% Offset used to map from ThermMech grid (positive down) to sea-level at 0 and positive up.
LEMpar.offset = sticky_layer;
LEMpar.sea_level = LEMpar.offset-water_lev;  % Define sea_level for surface processes and marine deposition, if water_lev is empty, it will be empty 
LEMPar.max_sea_level = LEMpar.sea_level;     % Define the absolut max of sea level, important for dynamic changes driven by max_depth.
LEMpar.max_depth = 500;  % Define maximum depth (m), this will cause the LEMpar.sea_level and water_lev to change accordingly 
                         % set to very high if no control is desired                        

% Apply surface processes model: 'none','line','fastscape'
% Make sure you edit the fastscape setup at the bottom of this file if
% using 'fastscape'
LEMpar.apply_surfacepro = 'none';

% Define erosional parameters
LEMpar.Kf = 0; % was 5e-6, bedrock river incision rate meters^(1-2m)/yr
LEMpar.Kd = 0;    % was 1, bedrock transport coefficient (diffusivity) meters^2/yr
LEMpar.Ks = 0;   % was 20, marine sediment transport coefficient
LEMpar.Fd = 0;    % was 2, Effeciency of fluvial sediment deposition (G in fastscape)
                  % You can turn off fluvial deposition by setting to zero.  
             
LEMpar.Ld = false ; % Deposit sediment in lakes true/false. 
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
% gridt(2,:) = sticky_layer;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% initial topography to hold the water
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% left and right side of the wedge topo
topo_x = gridt(1,:);
i_left = topo_x >= xl_1 & topo_x < wedge_center;
i_right = topo_x >= wedge_center & topo_x < xl_2;
i_flat = topo_x < xl_1| topo_x >= xl_2 ;

gridt(2, i_flat) =  sticky_layer;
gridt(2, i_left) =  sticky_layer  + tand_d*(gridt(1, i_left) - xl_1);
gridt(2, i_right) =  sticky_layer  + tand_d*(xl_2 - gridt(1, i_right) );

% Set-up eroisonal surrace and stratigraphy horizons

% Erosional surface, tracks surface as if no erosion/deposition
E = gridt(2,:);

% Age interval
dt_horizon = 200e3;  % yr
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

% Generate paleo topography
if strcmpi(LEMpar.apply_surfacepro,'line')
    paleo_topo = true;
    % Generate paleotopography palteau runing for a 0.2 Myrs
    nstep = ceil(200e3/LEMpar.dt_max);
    h = 1000 + zeros(size(gridt(2,:)));
    h(1) = 0;
    h(end) = 0;
    b =  h; 
    for i = 1:nstep
        [~,h,b,~] = ErosionDeposition1D(LEMpar.dt_max,LEMpar.dt_max,gridt(1,:),h,h,LEMpar.Ld,LEMpar.Fd,LEMpar.Kf,LEMpar.Kd);
    end
end
% ------------------------------------------------------------------------
%% FastScape Initialization and Setup
% -------------------------------------------------------------------------

if strcmpi(LEMpar.apply_surfacepro,'fastscape')

    % Setup the python environment for FastScape using pyenv
    pe = pyenv;
    if pe.Status == 'NotLoaded'
        pyenv('Version','/Users/rmoucha/.pyenv/shims/python')
        pyenv("ExecutionMode","OutOfProcess")
    end

    % Import libraries
    LEMpar.NP = py.importlib.import_module('numpy');
    LEMpar.FS = py.importlib.import_module('fastscapelib_fortran');

    % Intialize FastScape

    LEMpar.FS.fastscape_init();

    % Setup FastScape environment and grid
    LEMpar.xl = xsize0;
    LEMpar.yl = 100e3;                       % 100 km y-grid
    LEMpar.nx = tnum;
    LEMpar.ny = ceil(LEMpar.yl/tstp + 1);    % same dy as dx
    LEMpar.nn = LEMpar.nx*LEMpar.ny;
    
    LEMpar.FS.fastscape_set_nx_ny(LEMpar.nx,LEMpar.ny);
    LEMpar.FS.fastscape_setup();
    LEMpar.FS.fastscape_set_xl_yl(LEMpar.xl,LEMpar.yl);

    % Create grid
    zero_vector = zeros(1,LEMpar.nn);

    % Surfrace processes parameters
    Kf = LEMpar.Kf + zero_vector;
    Kd = LEMpar.Kd + zero_vector;
    Kf = LEMpar.NP.array(Kf);
    Kd = LEMpar.NP.array(Kd);
    LEMpar.FS.fastscape_set_erosional_parameters(Kf,LEMpar.kfsed,LEMpar.m,LEMpar.n,Kd,LEMpar.kdsed,LEMpar.Fd,LEMpar.Fd,LEMpar.expp);

    LEMpar.sea_level = -100; % Sea level in meters
    LEMpar.p1 = 0; % Sand porosity
    LEMpar.p2 = 0; % shale porosity
    LEMpar.z1 = 0; % e-folding depth for sand
    LEMpar.z2 = 0; % e-folding depth for shale
    LEMpar.r = 1; % sand-shale ratio at source
    LEMpar.L = 100; % averaging depth for solving shale-sand equation
    LEMpar.kds1 = 200; % marine transport coefficient for sand
    LEMpar.kds2 = 200; % marine transport coefficient for shale

    LEMpar.FS.fastscape_set_marine_parameters( LEMpar.sea_level, LEMpar.p1, LEMpar.p2, LEMpar.z1, LEMpar.z2, LEMpar.r, LEMpar.L, LEMpar.kds1, LEMpar.kds2)

    % Set the topography
    hinit = LEMpar.offset - gridt(2,:);
    hinit = hinit + 1000;
    hinit = repmat(hinit,LEMpar.ny,1);
    rng(1)
    hinit = hinit + randn(size(hinit));
    hinit(:,1) = 0;   % Set the base level
    hinit(:,end) = 0;
    hinit = reshape(hinit',1,LEMpar.nn);
    LEMpar.hfinal = LEMpar.NP.array(hinit);      % Containers allocated for later use
    LEMpar.bfinal = LEMpar.NP.array(hinit);      % Containters allocated for later use
    LEMpar.lake_depth = LEMpar.NP.array(zero_vector);  % Containters allocated for later use
    LEMpar.FS.fastscape_init_h(LEMpar.hfinal);

    % Precipitation
    p = LEMpar.p * ones(1,LEMpar.nn);
    p = LEMpar.NP.array(p);
    LEMpar.FS.fastscape_set_precip(p);

    % Set boundary conditions
    bc = 101; % Left and Right fixed, Top and Bottom cyclic
    LEMpar.FS.fastscape_set_bc(int32(bc));

    % Generate paleotopography palteau runing for a 0.2 Myrs
    nstep = ceil(200e3/LEMpar.dt_max);
    LEMpar.FS.fastscape_set_dt(LEMpar.dt_max);

    for i = 1:nstep
        % Run FastScape
        LEMpar.FS.fastscape_execute_step();
    end

    % Retrieve the topography and basement from FastScape
    LEMpar.FS.fastscape_copy_h(LEMpar.hfinal);
    LEMpar.FS.fastscape_copy_basement(LEMpar.bfinal);

    % Take the averages of topography and basement in the FastScape model
    % y-space, since we don't transpose, we take mean along dimension 2
    h = mean(reshape(double(LEMpar.hfinal),LEMpar.nx,LEMpar.ny),2);
    b = mean(reshape(double(LEMpar.bfinal),LEMpar.nx,LEMpar.ny),2);
    paleo_topo = true;

end

paleo_topo = false; 

if paleo_topo
    % Bring the topography and basement back into the ThermoMech model y-space
    gridt(2,:) = LEMpar.offset-h;
    H(1,:) = LEMpar.offset-b;

    for mm1 = 1:marknum

        % Erosion-sedimentation
        % Find topography node for the marker
        xn=double(int16((MX(mm1)-gridt(1,1))/tstp-0.5))+1;
        if (xn<1)
            xn=1;
        end
        if (xn>tnum-1)
            xn=tnum-1;
        end
        % Save horizontal index
        TXN(mm1)=xn;

    end

    % Compute relative distance to topography node
    dx = (MX-gridt(1,TXN)')/tstp;
    % Compute topograhy elevation above the marker
    dy = gridt(2,TXN)'.*(1-dx)+gridt(2,TXN+1)'.*dx;

    % water/air to rock transformation
    k = (MI == 1 | MI == 2) & MY > dy;
    MI(k) = 4;   % Change marker type to basalt = 4 or crust = 6
    MRAT(k) = 1; % Reset strain rate ratio
    MGII(k) = 0; % Reset strain

    % Rocks to air transformation
    k = MI > 1 & MI~=2 & MY < dy;

    MI(k) = 1;   % Change marker type
    MRAT(k) = 1; % Reset strain rate ratio
    MGII(k) = 0; % Reset strain

end

