% Define Material Properties, Build Geometry, Construct Boundary Conditions
%

%% Define Material Properties by loading a file

llnl_mat
matprop_copy = [outpath,'/matprop.m'];
copyfile('llnl_mat.m',matprop_copy);

% -------------------------------------------------------------------------
%% Defining Phase Transformations
% (1) Basalt to eclogite
% (2) 410 km
% (3) 660 km
% -------------------------------------------------------------------------
MPHASE = zeros(9,4);

% Eclogite (prograde and retrograde)

% Possible Phase Transformation types:
% false = no phase transformationsb
% 'depth' = pressure based according MPHASE(4,:)
% 'thermodynamic' = based on thermodynamic criteria (Gibbs Free Energy)
%                   with blueshicst cut-off
% 'kinetic' = based on kinetics

RxType =  'kinetic';
phase_function = 'univariate';   % Univeriate or hyperbollic

MPHASE(1,1:2) = [175e3 175e3];    % Kinetic Activation energy [prograde retrograde]
MPHASE(2,1:2) = [1e-6 1e-6];      % Kinetic Pre-exponential factor [prograde retrograde]
MPHASE(3,1:2) = 4;                % Avrami exponent 
                                  % = 4 means rate of transformation is not constant, but starts 
                                  %     at zero and ends at zero
                                  % = 1 reaction is fast at first than slow
MPHASE(4,1:2)= 1.5e9;             % Fixed Pressure phase transition (Pa)
MPHASE(5,1:2) = 200;             % Density Change (kg/m^3)

% 410 km
MPHASE(4,3) = 1.5;              % Clapeyron slope (MPa/K)
% MPHASE(5,3) = 206;              % Density change (kg/m^3)
MPHASE(5,3) = 0;              % Density change (kg/m^3)
MPHASE(6,3) = 13.5;             % Reference phase change pressure (GPa)
MPHASE(7,3) = 1810;             % Reference phase change temperature (K)
MPHASE(8,3) = 1;                % Viscosity increase factor
MPHASE(9,3) = 0.1;              % Phase transition range for hyperbolic (0.1 GPa = 20 km)

% 660 km
MPHASE(4,4) = -1.3;             % Clapeyron slope (MPa/K)
MPHASE(5,4) = 322;              % Density change (kg/m^3)
MPHASE(6,4) = 23;               % Reference phase change pressure (GPa)
MPHASE(7,4) = 1940;             % Reference phase change temperature (K)
MPHASE(8,4) = 1;                % Viscosity increase factor
MPHASE(9,4) = 0.1;              % Phase transition range for hyperbolic (0.1 GPa = 20 km)

% -------------------------------------------------------------------------
%% Defining lithological temperature model constants
% -------------------------------------------------------------------------

% For forced convergence

% Sticky Layer thickness water or air (m)
sticky_layer = 15000;

% Water Level Tracking (m) set to empty if no water
% Water depth = sticky_layer - water_lev
water_lev = 7000;

% % Sticky Layer thickness water or air (m)
% sticky_layer = 10000;
% 
% % Water Level Tracking (m) set to empty if no water
% % Water depth = sticky_layer - water_lev
% water_lev = 0;

% Sediment Thickness (m)
sed_thick = 1e3;

% Thermal diffusivity of the mantle, m^2/s
kappa = 1e-6;    

% Oceanic lithosphere thickness (1000oC isotherm) (m)
isoLAB = 1000-tabsolute;

% Age of Oceanic Plate (s)
global plate_age_old
plate_age_old = plate_age_old*yr2sec;  % Subducting

% Age of Oceanic Plate (s)
plate_age_young = 50e6*yr2sec; % Overriding

%------------------------------------------------------
% Setup subduction weak zone paramters 
%------------------------------------------------------

% Slab dipping angle (to the right)
dip_angle = 30;  % degrees

% Slope of subduction zone
m = tand(dip_angle);

% Diping Weak Zone, top of zone at bottom sediments, corresponds to slab

% Width of weak zone (normal to dip) 
weak_width = 5e3;                    % (m)
x1 = init_trenchpos;                           % Left Weak zone position (m)
x2 = x1 + weak_width/sind(dip_angle);  % Right Weak zone

% Intercept
b1 = sticky_layer + sed_thick - m*x1;
b2 = sticky_layer + sed_thick - m*x2;

%--------------------------------------------------------------------------
% Setup continent parameters
%--------------------------------------------------------------------------

% Continent-Shelf Transition (m)
cont_shelf = 50e3;

% Continent level (m), use this to set continent above or at sea level
cont_lev = water_lev-1000;    % <--- adjust for isostasy

% Upper Continental Crust thickness (m)
upper_cont_thick = 20000;

% Lower Continental Crust thickness (m)
lower_cont_thick = 15000;

% Mantle lithosphere thickness (m)
cont_mantle_lith_thick = 165000;

% Lithosphere thickness (m)
cont_litho_thick = cont_mantle_lith_thick + upper_cont_thick + lower_cont_thick;

% Bottom of upper crust continental
bot_ccrust1 = cont_lev + upper_cont_thick;

% Bottom of lower crust continental
bot_ccrust2 = bot_ccrust1 + lower_cont_thick;

% Continental lithosphere asthenosphere boundary (m) 
contLAB = cont_lev + cont_litho_thick;

%--------------------------------------------------------------------------
%Setup shelf parameters
%--------------------------------------------------------------------------

% Continent-Ocean Transition (m)
cont_ocean = cont_shelf + 1450e3;

% Shelf elevation at edge of shelf
shelf_lev = cont_lev +4e3;

% Upper Continental Crust thickness (m) at right edge of shelf
upper_cont_thick_right = 12e3;

% Lower Continental Crust thickness (m) at right edge of shelf
lower_cont_thick_right = 9e3;

% Mantle lithosphere thickness (m) at right edge of shelf
%cont_mantle_lith_thick_right = 79e3;
cont_mantle_lith_thick_right = 59e3;

% Lithosphere thickness (m) at right edge of shelf
cont_litho_thick_right = cont_mantle_lith_thick_right + upper_cont_thick_right + lower_cont_thick_right;

% Bottom of upper crust continental at right edge of shelf
bot_ccrust1_right = water_lev + 4e3 + upper_cont_thick_right;
%bot_ccrust1_right = sticky_layer + 2000 + upper_cont_thick_right;

% Bottom of lower crust continental at right edge of shelf
bot_ccrust2_right = bot_ccrust1_right + lower_cont_thick_right;

% Continental lithosphere asthenosphere boundary (m) at right edge of shelf
contLAB_right = shelf_lev + cont_litho_thick_right;

%--------------------------------------------------------------------------
%% Initial Temperature parameters
%--------------------------------------------------------------------------

% Convert deg-C to K
ttop = ttop - tabsolute;

%--------------------------------------------------------------------------
% for updated continental geotherm

% Temperatures at MOHO and LAB
tcontMoho = 550 - tabsolute;
tcontLAB = 1300 - tabsolute;

% % Geotherm parameters (K): slope and intercepts for crust and mantle
% % lithosphere
% m_crust = (tmoho - ttop) / (crust_thick);
% b_crust = tmoho - m_crust * crust_thick;
m_mlith = (tcontLAB - tcontMoho) / (contLAB - bot_ccrust2);
b_mlith = tcontMoho - m_mlith * bot_ccrust2;

% Calculate intersection depth of mantle lithosphere geotherm with adiabat
adiabat_depth = (b_mlith - (tpotential-tabsolute))/(tgrad - m_mlith);

% % Temperature at bottom of continental lithosphere
% tcontLAB = tpotential - tabsolute + tgrad*(contLAB - cont_lev);
%--------------------------------------------------------------------------

% Temperature at bottom of model space
mantle_thick = ysize-sticky_layer;
tbottom = tpotential - tabsolute + tgrad*mantle_thick;
thermLAB = 250e3;
% tbottom = 1800;

% Half-space cooling geotherm
z = linspace(0,mantle_thick,20000);
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

% Bottom of sediments ocean
bot_osed = sticky_layer + sed_thick;
bot_sed = sed_thick + sticky_layer;

% Upper Ocean Crust thickness (m)
upper_ocrust_thick = 3e3;
bot_ocrust1 = bot_osed + upper_ocrust_thick;

% lower Ocean Crust thickness (m)
lower_ocrust_thick = 4e3;
bot_ocrust2 = bot_ocrust1 + lower_ocrust_thick;

mantle_lith_thick_old = oceanic_lith_thick_old - (sed_thick + upper_ocrust_thick + lower_ocrust_thick);
mantle_lith_thick_young = oceanic_lith_thick_young - (sed_thick + upper_ocrust_thick + lower_ocrust_thick);

% Ocean lithosphere asthenosphere boundary (m) 
oceanLAB_old = sticky_layer + oceanic_lith_thick_old;
oceanLAB_young = sticky_layer + oceanic_lith_thick_young;

% -------------------------------------------------------------------------
%% Defining lithological structure of the model
% -------------------------------------------------------------------------

% Marker counter
mm1=0;
rng('default') % seed

% Initial trench position (m)
%TrenchX = 750e3;
TrenchX = init_trenchpos;
TrenchX_0 = TrenchX + 50e3;

% % Intial upper portion of slab tip depth, add sticky layer later (m)
% DepthY = sticky_layer + 200e3;
% 
% % Initial radius of curvature with respect to surface (m)
% Rinit = 500e3;
% 
% 
% % Initial radius of curvature of weak zone with respect to surface (m)
% weakz_thick = 20e3;
% Rweak = Rinit + weakz_thick;

% Loop over makers initial loop
for mm1 = 1:marknum
    
    % Dipping weak zone x-coordinates
    LX = (MY(mm1)-b1)/m;    %left bounded
    RX = (MY(mm1)-b2)/m;    %right bounded x
    
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
    if(MY(mm1)>sticky_layer && MY(mm1)<= bot_sed && MX(mm1)<=(xsize-500e3))
        MI(mm1)=3;
    end
   
    % Basaltic crust
    if(MY(mm1)>bot_osed && MY(mm1)<=bot_ocrust1 && MX(mm1)<=(xsize-500e3))
        MI(mm1)=4;
    end
    
    % Gabbroic crust
    if(MY(mm1)>bot_ocrust1 && MY(mm1)<=bot_ocrust2 && MX(mm1)<=(xsize-500e3))
        MI(mm1)=5;
    end
    
    %--------------------------------------------------------
    % Subducting Ocean Plate
    %--------------------------------------------------------
    
    % Mantle lithosphere
    if(MY(mm1)>bot_ocrust2 && MY(mm1)<=oceanLAB_old && MX(mm1)<=LX)
        MI(mm1) = 8;
    end 
    
    %--------------------------------------------------------
    % Overriding Ocean Plate
    %--------------------------------------------------------
    
    % Mantle lithospre
    if(MY(mm1)>bot_ocrust2 && MY(mm1)<=oceanLAB_young && MX(mm1)>=RX && MX(mm1)<=(xsize-500e3))
        MI(mm1) = 8;
    end
 
    

     global microcontinent_right
    %--------------------------------------------------------
    % Microcontinent
    %-------------------------------------------------------- 
     if(MY(mm1)>shelf_lev && MY(mm1)<=bot_ccrust1_right && MX(mm1)>=RX &&  MX(mm1)<=(xsize-microcontinent_right))
        MI(mm1)=6;
     end
     
     if(MY(mm1)>bot_ccrust1_right && MY(mm1)<=bot_ccrust2_right && MX(mm1)>=RX &&  MX(mm1)<=(xsize-microcontinent_right))
        MI(mm1)=7;
     end
     
     if(MY(mm1)>bot_ccrust2_right && MY(mm1)<=contLAB_right && MX(mm1)>=RX &&  MX(mm1)<=(xsize-microcontinent_right))
        MI(mm1)=9;
     end
     
%     %--------------------------------------------------------
%     % Fluid weakened part of overriding plate
%     %--------------------------------------------------------
%      if(MY(mm1)>bot_ccrust2_right && MY(mm1)<=contLAB_right && MX(mm1)>=(4350e3) &&  MX(mm1)<=(4400e3))
%         MI(mm1)=11;
%      end
    
    
    %--------------------------------------------------------
    % Weak Zone
    %--------------------------------------------------------
    
     if(MY(mm1)>bot_osed && MY(mm1)<=oceanLAB_old && MX(mm1)>=LX && MX(mm1)<=RX)
        MI(mm1)=11;
     end

end

% %--------------------------------------------------------------------------
% % Material Properties of subducting portion of the slab
% %--------------------------------------------------------------------------
%         
% % Arc angle
% beta = asind((Rinit-DepthY)./Rinit); 
% theta = degtorad(linspace(beta,90,100));
% 
% % Asemble the top arc
% [x1,y1] = pol2cart(theta,Rinit);
% 
% % Set coordinate to be positive with depth
% y1 = abs(y1-Rinit);
% 
% %--------------------------------------------------------------------------
% % Mantle Lithosphere
% %--------------------------------------------------------------------------
%  
% % Asemble the bottom arc
% [x2,y2] = pol2cart(theta,Rinit-oceanic_lith_thick_old);
% 
% % Set coordinate to be positive with depth
% y2 = abs(y2-Rinit);
% 
% % Assemble the polygon
% xv = [x1 fliplr(x2) x1(1)] + TrenchX;
% yv = [y1 fliplr(y2) y1(1)] + sticky_layer;
% 
% % Search for markers in polygon
% [in, on] = inpolygon(MX,MY,xv,yv);
% in = in | on;
% 
% % All oceanic lithosphere
% MI(in) = 8;
% 
% %--------------------------------------------------------------------------
% % Lower Crust
% %--------------------------------------------------------------------------
%  
% % Asemble the bottom arc for Gabroic lower oceanic crust
% [x2,y2] = pol2cart(theta,Rinit-lower_ocrust_thick-upper_ocrust_thick-sed_thick);
% 
% % Set coordinate to be positive with depth
% y2 = abs(y2-Rinit);
% 
% % Assemble the polygon
% xv = [x1 fliplr(x2) x1(1)] + TrenchX;
% yv = [y1 fliplr(y2) y1(1)] + sticky_layer;
% 
% % Search for markers in polygon
% [in, on] = inpolygon(MX,MY,xv,yv);
% in = in | on;
% 
% % Lower crust
% MI(in) = 5;
% 
% %--------------------------------------------------------------------------
% % Upper Crust
% %--------------------------------------------------------------------------
% 
% % Asemble the bottom arc for basalic upper oceanic crust
% [x2,y2] = pol2cart(theta,Rinit-upper_ocrust_thick-sed_thick);
% 
% % Set coordinate to be positive with depth
% y2 = abs(y2-Rinit);
% 
% % Assemble the polygon
% xv = [x1 fliplr(x2) x1(1)] + TrenchX;
% yv = [y1 fliplr(y2) y1(1)] + sticky_layer;
% 
% % Search for markers in polygon
% [in, on] = inpolygon(MX,MY,xv,yv);
% in = in | on;
% 
% % Upper Crust
% MI(in) = 4;
% 
% %--------------------------------------------------------------------------
% % Sediments
% %--------------------------------------------------------------------------
% 
% % Asemble the bottom arc for sediment (weak zone)
% [x2,y2] = pol2cart(theta,Rinit-sed_thick);
% 
% % Set coordinate to be positive with depth
% y2 = abs(y2-Rinit);
% 
% % Assemble the polygon
% xv = [x1 fliplr(x2) x1(1)] + TrenchX;
% yv = [y1 fliplr(y2) y1(1)] + sticky_layer;
% 
% % Search for markers in polygon
% [in, on] = inpolygon(MX,MY,xv,yv);
% in = in | on;
% 
% % All oceanic lithosphere
% MI(in) = 3;
% 
% %--------------------------------------------------------------------------
% % Weak zone
% %--------------------------------------------------------------------------
% 
% % Arc angle
% betaw = asind((Rinit-contLAB_right+sticky_layer)./Rweak); % Bottom angle
% thetaw = degtorad(linspace(betaw,90,100));
% 
% % Asemble the top arc
% [x1w,y1w] = pol2cart(thetaw,Rweak);
% 
% % Set coordinate to be positive with depth
% y1w = Rinit-y1w;
% 
% % Asemble the bottom arc for sediment (weak zone)
% [x2,y2] = pol2cart(thetaw,Rinit);
% 
% % Set coordinate to be positive with depth
% y2 = Rinit-y2;
% 
% % Assemble the polygon
% xv = [x1w fliplr(x2) x1w(1)] + TrenchX;
% yv = [y1w fliplr(y2) y1w(1)] + sticky_layer;
% 
% % Search for markers in polygon bounded by bottom of oceanic crust and the
% % LAB of young lithosphere
% [in, on] = inpolygon(MX,MY,xv,yv);
% 
% % Add extra sediments to simulate trench
% in1 = (in | on) & MY >= bot_sed & MY < bot_ocrust1;
% MI(in1) = 3;
% 
% in2 = (in | on) & MY >= bot_ocrust1 & MY  <= contLAB_right;
% % All oceanic lithosphere and lower curst changed to hydrated mantle
% MI(in2) = 11;

%--------------------------------------------------------------------------
% Initial temperature structure
%--------------------------------------------------------------------------

for mm1 = 1:marknum   
    
    
    % Dipping weak zone x-coordinates
    LX = (MY(mm1)-b1)/m;    %left bounded
    RX = (MY(mm1)-b2)/m;    %right bounded x
    
    % Adiabatic Temperature Gradient in the Mantle
    MTK(mm1)=tbottom-tgrad*(ysize-MY(mm1));
    
    % Sticky air or sticky water
    if(MI(mm1) == 1 || MI(mm1) == 2)
        MTK(mm1) = ttop;
    end


    
    %--------------------------------------------------------------------------

    % Overriding (young) plate
    if ( MY(mm1) > sticky_layer && MY(mm1) <= thermLAB && MX(mm1) > RX && MX(mm1) <= (xsize-500e3))
        z = MY(mm1) - sticky_layer;
        T_adiabat = tpotential - tabsolute + tgrad*z;
        MTK(mm1) = ttop + (T_adiabat-ttop).*erf(z./(2*sqrt(kappa*plate_age_young)));
    end
    
    % Downgoing (old) plate
    if ( MY(mm1) > sticky_layer && MY(mm1) <= thermLAB && MX(mm1) <= RX)
        z = MY(mm1) - sticky_layer;
        T_adiabat = tpotential - tabsolute + tgrad*z;
        MTK(mm1) = ttop + (T_adiabat-ttop).*erf(z./(2*sqrt(kappa*plate_age_old)));
    end
    
    %--------------------------------------------------------------------------
    % Microcontinent temperature
    %--------------------------------------------------------------------------
    
    % % Geotherm parameters (K): slope and intercepts for crust and mantle
    m_mlith_right = (tcontLAB - tcontMoho) / (contLAB_right - bot_ccrust2_right);
    b_mlith_right = tcontMoho - m_mlith_right * bot_ccrust2_right;

    % Calculate intersection depth of mantle lithosphere geotherm with adiabat
    adiabat_right = (b_mlith_right - (tpotential-tabsolute))/(tgrad - m_mlith_right); 
    
    
    % Continental geotherm is linear from surface to Moho
    if( MY(mm1) > shelf_lev && MY(mm1) <= bot_ccrust2_right && MX(mm1)>=RX &&  MX(mm1)<=(xsize-microcontinent_right))
        MTK(mm1) = ttop + (tcontMoho-ttop)*(MY(mm1)-shelf_lev)/(bot_ccrust2_right-shelf_lev);
        %MTK(mm1) = 550;
    end
    
    % Continental geotherm is linear from Moho to LAB (adiabat depth)
    if( MY(mm1) > bot_ccrust2_right && MY(mm1) <= adiabat_right && MX(mm1)>=RX &&  MX(mm1)<=(xsize-microcontinent_right))
        MTK(mm1) = tcontMoho + (tcontLAB-tcontMoho)*(MY(mm1)-bot_ccrust2_right)/(contLAB_right-bot_ccrust2_right);
        %MTK(mm1) = 1300;
    end
    
    
end

%--------------------------------------------------------------------------
% Plume parameters
%--------------------------------------------------------------------------
add_plume = false;

Texcess = 300;

% Plume head size and position (m)
rplume = 300e3;
xplume = 2500e3;
yplume = 1100e3;

% %--------------------------------------------------------------------------
% % Thermal Structure of subducting portion of the slab
% %--------------------------------------------------------------------------
% 
% % Asemble the bottom arc (remove sticky_layer)
% [x2,y2] = pol2cart(theta,Rinit-thermLAB);
% 
% % Set coordinate to be positive with depth
% y2 = abs(y2-Rinit);
% 
% % Assemble the polygon
% xv = [x1 fliplr(x2) x1(1)] + TrenchX;
% yv = [y1 fliplr(y2) y1(1)] + sticky_layer;
% 
% % Search for markers in polygon
% [in, on] = inpolygon(MX,MY,xv,yv);
% in = in | on;
% 
% z = Rinit - sqrt((MX(in)-TrenchX).^2 + (MY(in)-Rinit-sticky_layer).^2);
% Tz = tpotential - tabsolute + tgrad*z;
% MTK(in) = ttop + (Tz - ttop).*erf(z/(2*sqrt(kappa*plate_age_old)));

% Save Number of markers marknum=mm1
marknum = mm1;



% ------------------------------------------------------------------------
%% Topography and Surface Processes
% -------------------------------------------------------------------------

% Use a surface process timestep <5000 yrs.
LEMpar.dt_max = 2000; %yrs

% Apply surface processes model true/false
LEMpar.apply_surfacepro = false;

global Kf 
if isempty(Kf)
    Kf = 1e-5;
end
global Kd 
if isempty(Kd)
    Kf = 1;
end

% LEM parameters
LEMpar.Kf = 5e-6; % bedrock river incision rate meters^(1-2m)/yr
LEMpar.Kd = 1;    % bedrock transport coefficient (diffusivity) meters^2/yr
LEMpar.Ks = 500;  % marine sediment transport coefficient
LEMpar.Fd = 1;    % Effeciency of fluvial sediment deposition. 
                  % You can turn off fluvial deposition by setting to zero.  
LEMpar.sea_level = water_lev;  % Define sea_level for surface processes and marine deposition  
LEMpar.Lcrit = 5e3; % Marine sediment deposition critical distance from shore (m)
LEMpar.Ld = true;  % Deposit sediment in lakes true/false. 


% Desired topography model resolution
global tstp
if isempty(tstp)
    tstp = 250; % (m)
end

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


% -------------------------------------------------------------------------
%% Initial elevation for topography profile
% -------------------------------------------------------------------------

% Intial elevation at bottom of sticky layer
% gridt(2,:) = sticky_layer;

% % Set initial continent elevation
% k = gridt(1,:) >= ix & gridt(1,:);
% gridt(2,k) = cont_lev;

% % Set initial continent elevation
% k = gridt(1,:) <= cont_shelf;
% gridt(2,k) = cont_lev;
% 
% k = gridt(1,:) <= cont_ocean & gridt(1,:) > cont_shelf;
% gridt(2,k) =(m_ccrust_top* gridt(1,k) +b_ccrust_top);

% k = gridt(2,:) > sticky_layer;
% gridt(2,k) = sticky_layer;

% k = gridt(1,:) <= cont_ocean;
% gridt(2,k) = cont_lev;

% Save initial topography
% topotime(1,1) = 0;
% topohigh(1,:) = gridt(2,:);
% if ~isempty(water_lev)
%     topo_lev=zeros(ntimestep+1,tnum);
%     topowater(1,1) = water_lev;
% end