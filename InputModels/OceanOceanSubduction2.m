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

% Thermal diffusivity of the mantle, m^2/s
kappa = 1e-6;    

% Oceanic lithosphere thickness (1000oC isotherm) (m)
isoLAB = 1100-tabsolute;

% Bottom of thermal lithosphere
ThermLith = 1300-tabsolute;

% Age of Oceanic Plate (s)
plate_age_old = 60e6*yr2sec;  % Subducting

% Age of Oceanic Plate (s)
plate_age_young = 40e6*yr2sec; % Overriding

% 

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
T_old = ttop + (T_adiabat-ttop).*erf(z./(2*sqrt(kappa*plate_age_old)));
T_young = ttop + (T_adiabat-ttop).*erf(z./(2*sqrt(kappa*plate_age_young)));

% Seimic oceanic lithosphere thickness at isoLAB
k = find(T_old>=isoLAB,1);
if isempty(k)
    error('Model space exceeds isoLAB')
else
    oceanic_lith_thick_old = z(k);
end
k = find(T_young>=isoLAB,1);
if isempty(k)
    error('Model space exceeds isoLAB')
else
    oceanic_lith_thick_young = z(k);
end

% Thermal oceanic lithosphere depth at ThermLith
k = find(T_old>ThermLith,1);
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
oceanLAB_old = sticky_layer + oceanic_lith_thick_old;
oceanLAB_young = sticky_layer + oceanic_lith_thick_young;
oceanic_lith_base = sticky_layer + oceanic_lith_base;

% -------------------------------------------------------------------------
%% Defining lithological structure of the model
% -------------------------------------------------------------------------

% Marker counter
mm1=0;
rng('default') % seed

% Initial trench position (m) -- middle of model space
TrenchX = 1500e3;%xsize0/4;

% Intial upper portion of slab tip depth, add sticky layer later (m)
DepthY = oceanic_lith_thick_old+30e3;

% Initial radius of curvature with respect to surface (m)
Rinit = 500e3;

% Initial radius of curvature of weak zone with respect to surface (m)
weakz_thick = 10e3;
Rweak = Rinit + weakz_thick;

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
    if(MY(mm1)>bot_ocrust2 && MY(mm1)<=oceanLAB_old && MX(mm1)<=TrenchX)
        MI(mm1) = 8;
    end
    
    %--------------------------------------------------------
    % Overriding Young Ocean Plate
    %--------------------------------------------------------
        
    % Mantle lithospre
    if(MY(mm1)>bot_ocrust2 && MY(mm1)<=oceanLAB_young && MX(mm1)>=TrenchX)
        MI(mm1) = 8;
    end
 
    % Weak zone
      if (MY(mm1)> bot_ocrust2 && MY(mm1)<=oceanLAB_young && MX(mm1)>= xsize0 - oceanic_lith_thick_young/2)
           MI(mm1) = 11;
      end
      
     if (MY(mm1)> bot_ocrust2 && MY(mm1)<=oceanLAB_old && MX(mm1)<=oceanic_lith_thick_old)
          MI(mm1) = 11;
     end
end

%--------------------------------------------------------------------------
% Material Properties of subducting portion of the slab
%--------------------------------------------------------------------------
        
% Arc angle
beta = asind((Rinit-DepthY)./Rinit); 
theta = deg2rad(linspace(beta,90,100));

% Asemble the top arc
[x1,y1] = pol2cart(theta,Rinit);

% Set coordinate to be positive with depth
y1 = Rinit-y1;

%--------------------------------------------------------------------------
% Mantle Lithosphere
%--------------------------------------------------------------------------
 
% Asemble the bottom arc
[x2,y2] = pol2cart(theta,Rinit-oceanic_lith_thick_old);

% Set coordinate to be positive with depth
y2 = Rinit-y2;

% Assemble the polygon
xv = [x1 fliplr(x2) x1(1)] + TrenchX;
yv = [y1 fliplr(y2) y1(1)] + sticky_layer;

% Search for markers in polygon
[in, on] = inpolygon(MX,MY,xv,yv);
in = in | on;

% All oceanic lithosphere
MI(in) = 8;

%--------------------------------------------------------------------------
% Lower Crust
%--------------------------------------------------------------------------
 
% Asemble the bottom arc for Gabroic lower oceanic crust
[x2,y2] = pol2cart(theta,Rinit-lower_ocrust_thick-upper_ocrust_thick-sed_thick);

% Set coordinate to be positive with depth
y2 = Rinit-y2;

% Assemble the polygon
xv = [x1 fliplr(x2) x1(1)] + TrenchX;
yv = [y1 fliplr(y2) y1(1)] + sticky_layer;

% Search for markers in polygon
[in, on] = inpolygon(MX,MY,xv,yv);
in = in | on;

% Lower crust
MI(in) = 5;

%--------------------------------------------------------------------------
% Upper Crust
%--------------------------------------------------------------------------

% Asemble the bottom arc for basalic upper oceanic crust
[x2,y2] = pol2cart(theta,Rinit-upper_ocrust_thick-sed_thick);

% Set coordinate to be positive with depth
y2 = Rinit-y2;

% Assemble the polygon
xv = [x1 fliplr(x2) x1(1)] + TrenchX;
yv = [y1 fliplr(y2) y1(1)] + sticky_layer;

% Search for markers in polygon
[in, on] = inpolygon(MX,MY,xv,yv);
in = in | on;

% Upper Crust
MI(in) = 4;

%--------------------------------------------------------------------------
% Sediments
%--------------------------------------------------------------------------

% Asemble the bottom arc for sediment (weak zone)
[x2,y2] = pol2cart(theta,Rinit-sed_thick);

% Set coordinate to be positive with depth
y2 = Rinit-y2;

% Assemble the polygon
xv = [x1 fliplr(x2) x1(1)] + TrenchX;
yv = [y1 fliplr(y2) y1(1)] + sticky_layer;

% Search for markers in polygon
[in, on] = inpolygon(MX,MY,xv,yv);
in = in | on;

% All oceanic lithosphere
MI(in) = 3;

%--------------------------------------------------------------------------
% Weak zone
%--------------------------------------------------------------------------

% Arc angle
betaw = asind((Rinit-oceanLAB_young+sticky_layer)./Rweak); % Bottom angle
thetaw = deg2rad(linspace(betaw,90,100));

% Asemble the top arc
[x1w,y1w] = pol2cart(thetaw,Rweak);

% Set coordinate to be positive with depth
y1w = Rinit-y1w;

% Asemble the bottom arc for sediment (weak zone)
[x2,y2] = pol2cart(thetaw,Rinit);

% Set coordinate to be positive with depth
y2 = Rinit-y2;

% Assemble the polygon
xv = [x1w fliplr(x2) x1w(1)] + TrenchX;
yv = [y1w fliplr(y2) y1w(1)] + sticky_layer;

% Search for markers in polygon bounded by bottom of oceanic crust and the
% LAB of young lithosphere
[in, on] = inpolygon(MX,MY,xv,yv);

% Add extra sediments to simulate trench
in1 = (in | on) & MY >= bot_sed & MY < bot_ocrust1;
MI(in1) = 3;

in2 = (in | on) & MY >= bot_ocrust1 & MY  <= oceanLAB_young;
% All oceanic lithosphere and lower curst changed to hydrated mantle
MI(in2) = 11;

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
    if( MY(mm1) > sticky_layer && MY(mm1)<=oceanic_lith_base && MX(mm1) <= TrenchX)
       z = MY(mm1) - sticky_layer;
       MTK(mm1) = ttop + (MTK(mm1)-ttop).*erf(z/(2*sqrt(kappa*plate_age_old)));
    end
    
    % Oceanic geotherm for overriding plate of specified age
    % after Turcotte & Schubert (2002)
    if( MY(mm1) > sticky_layer && MY(mm1)<=oceanic_lith_base && MX(mm1) > TrenchX)   
        z = MY(mm1) - sticky_layer;
        MTK(mm1) = ttop + (MTK(mm1)-ttop).*erf(z/(2*sqrt(kappa*plate_age_young)));
    end
    
end

%--------------------------------------------------------------------------
% Thermal Structure of subducting portion of the slab
%--------------------------------------------------------------------------

% Asemble the bottom arc (remove sticky_layer)
[x2,y2] = pol2cart(theta,Rinit-oceanic_lith_base+sticky_layer);

% Set coordinate to be positive with depth
y2 = Rinit-y2;

% Assemble the polygon
xv = [x1 fliplr(x2) x1(1)] + TrenchX;
yv = [y1 fliplr(y2) y1(1)] + sticky_layer;

% Search for markers in polygon
[in, on] = inpolygon(MX,MY,xv,yv);
in = in | on;

z = Rinit - sqrt((MX(in)-TrenchX).^2 + (MY(in)-Rinit-sticky_layer).^2);
Tz = tpotential - tabsolute + tgrad*z;
MTK(in) = ttop + (Tz - ttop).*erf(z/(2*sqrt(kappa*plate_age_old)));

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