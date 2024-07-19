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

MouchaMaterials_FaultCrust

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
mantle_lith = 120000;

% Lithosphere (LAB) (m) % Specifically, the Thermal lithosphere depth
LAB_depth = mantle_lith + crust_thick;

% Add a plume anomaly at given age?
% Add a plume anomaly at given age?
add_plume = 0.1;  % Ma,  set to false

% Plume head size and position (m)
plume_type = 'spherical'; %'gaussian'
rplumex = 150e3;
rplumey = 150e3;
xplume = xsize/2 + 75e3;%-175e3; % -20 km for 10 cm/yr %
yplume = 600e3;  % Will be set to ysize if exceeds ysize
%--------------------------------------------------------------------------
%% Initial Temperature parameters
%-------------------------------------------------------------------------

% Plume Excess Temperature [K]
Texcess = 300;

% Thermal diffusivity of the mantle, m^2/s
kappa = 1e-6;    

% Increase default mantle potential temperature oC
dtpotential = 0;%100;

% Convert deg-C to K
ttop = ttop - tabsolute;
tpotential = tpotential + dtpotential - tabsolute;

% Temperature at bottom of model space
tbottom = tpotential + tgrad*(ysize-sticky_layer);

% -------------------------------------------------------------------------
%% Defining lithological structure of the model
% -------------------------------------------------------------------------
iframe = 64;
iframe = num2str(iframe);
inpathold = '../AGU2018/Initial/vext30mmyr/Lith200dTp0_lrg/';

shift_x = 300e3;%650e3;

% LAB depth with sticky_layer is defined as the depth of 1200oC isotherm
% k = find(SE.T < 1200-tabsolute,1,'last');
% LAB_depth = SE.z(k) + sticky_layer;
LAB_depth = LAB_depth + sticky_layer;

% Use structure from previous run as initial input
load([inpathold,'markers_',iframe,'.mat'],'MX','MY','MI','MTK')

MX = MX - shift_x;

k = MX < gridx(1) | MX > gridx(xnum) | MY < gridy(1) | MY > gridy(ynum);

MX(k) = [];
MY(k) = [];
MI(k) = [];
MTK(k) = [];

marknum = length(MI);

% Initialize the rest
MXN=zeros(marknum,1);   % Horizontal index
MYN=zeros(marknum,1);   % Vertical index
MCXN=zeros(marknum,1);  % Horizontal central index
MCYN=zeros(marknum,1);  % Vertical central index

% -------------------------------------------------------------------------
%% Fill empty cells
% -------------------------------------------------------------------------

% Computing grid steps for basic nodes
xstp1 = gridx(2:xnum) - gridx(1:xnum-1);
ystp1 = gridy(2:ynum) - gridy(1:ynum-1);

% Computing grids and grid steps for Vx, Vy nodes
% Horizontal (for Vy)
gridcx=zeros(xnum+1,1);
xstpc1=zeros(xnum,1);
% Vertical (for Vx)
gridcy=zeros(ynum+1,1);
ystpc1=zeros(ynum,1);
% First and last nodes and steps (for external nodes)
% Horizontal (for Vy)
gridcx(1)=gridx(1)-xstp1(1)/2;
xstpc1(1)=xstp1(1);
gridcx(xnum+1)=gridx(xnum)+xstp1(xnum-1)/2;
xstpc1(xnum)=xstp1(xnum-1);
% Vertical (for Vx)
gridcy(1)=gridy(1)-ystp1(1)/2;
ystpc1(1)=ystp1(1);
gridcy(ynum+1)=gridy(ynum)+ystp1(ynum-1)/2;
ystpc1(ynum)=ystp1(ynum-1);

% Internal nodes
gridcx(2:xnum) = (gridx(2:xnum) + gridx(1:xnum-1))/2;
gridcy(2:ynum) = (gridy(2:ynum) + gridy(1:ynum-1))/2;

% Internal grid steps
xstpc1(2:xnum-1) = (gridx(3:xnum) - gridx(1:xnum-2))/2;
ystpc1(2:ynum-1) = (gridy(3:ynum) - gridy(1:ynum-2))/2;

for mm1 = 1:marknum
    
    %  xn    rho(xn,yn)--------------------rho(xn+1,yn)
    %           ?           ^                  ?
    %           ?           ?                  ?
    %           ?          dy                  ?
    %           ?           ?                  ?
    %           ?           v                  ?
    %           ?<----dx--->o Mrho(xm,ym)       ?
    %           ?                              ?
    %           ?                              ?
    %  xn+1  rho(xn,yn+1)-------------------rho(xn+1,yn+1)
    %
    % Define indexes for upper left node in the cell where the marker is
    % using bisection
    % Find horizontal index
    xnmin=1;
    xnmax=xnum;
    while ((xnmax-xnmin)>1)
        % !!! SUBTRACT 0.5 since int16(0.5)=1
        xn = round((xnmax+xnmin)/2-0.5);
        %xn=double(int16((xnmax+xnmin)./2-0.5));
        if(gridx(xn)>MX(mm1))
            xnmax=xn;
        else
            xnmin=xn;
        end
    end
    xn=xnmin;
    % Check horizontal index
    if (xn<1)
        xn=1;
    end
    if (xn>xnum-1)
        xn=xnum-1;
    end
    % Save horizontal index
    MXN(mm1)=xn;
    
    % Now define the index for the central nodes (Presure, SXX, etc..)
    if (MX(mm1) < gridcx(xn+1))
        xn = xn-1;
    end
    if (xn<1)
        xn=1;
    end
    if (xn>xnum-2)
        xn=xnum-2;
    end
    MCXN(mm1) = xn;
    
    % Find vertical index
    ynmin=1;
    ynmax=ynum;
    while ((ynmax-ynmin)>1)
        % !!! SUBTRACT 0.5 since int16(0.5)=1
        yn = round((ynmax+ynmin)./2-0.5);
        %yn=double(int16((ynmax+ynmin)./2-0.5));
        if(gridy(yn)>MY(mm1))
            ynmax=yn;
        else
            ynmin=yn;
        end
    end
    yn=ynmin;
    % Check vertical index
    if (yn<1)
        yn=1;
    end
    if (yn>ynum-1)
        yn=ynum-1;
    end
    % Save Vertical index
    MYN(mm1)=yn;
    
    % Now define the index for the central nodes (Presure, SXX, etc..)
    if (MY(mm1) < gridcy(yn+1))
        yn = yn-1;
    end
    if(yn<1)
        yn = 1;
    end
    if(yn>ynum-2)
        yn = ynum-2;
    end
    MCYN(mm1) = yn;
    
end

mx = []; my = []; mi = []; mxn = []; myn = []; newmarkers = 0;

% Because of the sticky-layer markers drift downwards over time, we must
% fill the empty half of the first row of cells. We do this by copying
% markes in the first row of the grid and moving them up by half the grid
% spacing
k = MYN == 1;
dy = gridy(2)/2;  % copy up by half grid
if sum(k)==0 % We have a problem, all the first row cells are empty, so we populat them with the next row
   k = MYN == 2;
   dy = gridy(2);  % copy by a full grid
end

if sum(k)==0
    error('First two rows are empty')
end

mx0 = MX(k); my0 = MY(k) - dy; mi0 = MI(k);
newmarkers0 = sum(k);

% Now we must fill any-empty cells with markers, which will most likely be
% asthenospheric (bottom left and right corners of the grid).
newmarkertype = 10; % Astenosphere marker
[mx,my,mxn,myn,mcxn,mcyn,newmarkers] = fillemptycell(gridx,gridy,gridcx,gridcy,MXN,MYN,mxcell,mycell,false);

zmarks = zeros(size(mx));
mi = newmarkertype+zmarks;
mtk = zmarks;

old = load([inpathold,'grids_',iframe,'.mat'],'G');

[X,Y]=meshgrid(old.G.gridx,old.G.gridy);
mtk0 = interp2(X,Y,old.G.tk1,mx0,my0);
mtk = interp2(X,Y,old.G.tk1,mx,my);

MX = [MX; mx0; mx];
MY = [MY; my0; my];
MTK = [MTK; mtk0; mtk];
MI = [MI; mi0; mi];

marknum = length(MI);

% Initialize the rest of the markers
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
plastyn = zeros(marknum,1); % plastic yielding occured
MPH410 = zeros(marknum,2);
MPH660 = MPH410;
MECL=zeros(marknum,5);   % Eclogite phase transformations

% -------------------------------------------------------------------------
%% Initial elevation for topography profile
% -------------------------------------------------------------------------

% Load topography from previous run as intial input
old = load([inpathold,'grids_',iframe,'.mat'],'gridt');
old.gridt(1,:) = old.gridt(1,:) - shift_x;

% Interpolate old topography to new gridt
gridt(2,:) = interp1(old.gridt(1,:),old.gridt(2,:),gridt(1,:));
k = isnan(gridt(2,:));
gridt(2,k) = sticky_layer;

