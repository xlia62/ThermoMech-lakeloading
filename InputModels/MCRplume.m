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
% Standard 14, Addons must start with 15 below
MatProp{1} = 'Sticky Air';
MatProp{2} = 'Sticky Water';
MatProp{3} = 'Sediments';
MatProp{4} = 'Upper Oceanic Crust';
MatProp{5} = 'Lower Oceanic Crust';
MatProp{6} = 'Upper Continental Crust';
MatProp{7} = 'Lower Continental Crust';
MatProp{8} = 'Oceanic Lithospheric Mantle';
MatProp{9} = 'Continental Lithospheric Mantle';
MatProp{10} = 'Asthenospheric Mantle';
MatProp{11} = 'Hydrated Mantle';
MatProp{12} = 'Upper Oceanic Plateau';
MatProp{13} = 'Lower Oceanic Plateau';
MatProp{14} = 'Plume Material';
% Add-ons HERE!

% MRHO  = density (kg/m3): RHO*[1-ALP*(T-273)]*[1+BET*(P-1e+5)], includes
%         melt density
% MFLOW = power-law: EPSILONii=AD*SIGMAii^n*exp[-(Ea+Va*P)/RT)
% MMU = shear modulus (Pa)
% MPL = Brittle/plastic strength (Pa): SIGMAyeild=C+sin(FI)*P
%         C=C0, FI=FI0 for strain<=GAM0
%         C=C0+(C1-C0)/(GAM1-GAM0)*(strain-GAM0), FI=FI0+(FI1-FI0)/(GAM1-GAM0)*(strain-GAM0) for GAM0<strain<GAM1
%         C=C1, FI=FI1 for strain>=GAM1
% MCP = heat capacity (J/K/kg)
% MKT = thermal conductivity (W/m/K): k=k0+a/(T+77)
% MHR = radiogenic heat production (W/m^3)
% MHE = shear and adiobatic heating efficiency:  take on values 0 or 1

% -------------------------------------------------------------------------
% 1 = Weak Layer ("sticky air")
% -------------------------------------------------------------------------
MRHO(1,1)=1.22;             % standard density, kg/m^3
MRHO(1,2)=3e-5;             % thermal expansion, 1/K
MRHO(1,3)=1e-11;            % compressibility, 1/Pa
MFLOW(1,1)=0;               % 0=constant viscosity
MFLOW(1,2)=1e+18;           % viscosity, Pa s
MMU(1,1)=1e+20;               % shear modulus, Pa
MPL(1,1)=0;                 % C0, Pa
MPL(1,2)=0;                 % C1, Pa
MPL(1,3)=0;                 % sin(FI0)
MPL(1,4)=0;                 % sin(FI1)
MPL(1,5)=0;                 % GAM0
MPL(1,6)=1;                 % GAM1
MCP(1,1)=3000;                % Cp, J/kg
MKT(1,1)=300;               % k0, W/m/K
MKT(1,2)=0;                 % a, W/m
MHR(1,1)=0;                 % radiogenic heat production, W/m^3
MHE(1,1)=0;                 % shear and adiobatic heating efficiency

% -------------------------------------------------------------------------
% 2 = Weak Layer ("sticky water")
% -------------------------------------------------------------------------
MRHO(2,1)=1000;             % standard density, kg/m^3
MRHO(2,2)=3e-5;             % thermal expansion, 1/K
MRHO(2,3)=1e-11;            % compressibility, 1/Pa
MFLOW(2,1)=0;               % 0=constant viscosity
MFLOW(2,2)=1e+18;           % viscosity, Pa s
MMU(2,1)=1e+20;               % shear modulus, Pa
MPL(2,1)=0;                 % C0, Pa
MPL(2,2)=0;                 % C1, Pa
MPL(2,3)=0;                 % sin(FI0)
MPL(2,4)=0;                 % sin(FI1)
MPL(2,5)=0;                 % GAM0
MPL(2,6)=1;                 % GAM1
MCP(2,1)=3000;                % Cp, J/kg
MKT(2,1)=300;               % k0, W/m/K
MKT(2,2)=0;                 % a, W/m
MHR(2,1)=0;                   % radiogenic heat production, W/m^3
MHE(2,1)=0;                 % shear and adiobatic heating efficiency

% -------------------------------------------------------------------------
% 3 = Sediments
% -------------------------------------------------------------------------
MRHO(3,1)=2600;             % standard density, kg/m^3
MRHO(3,2)=3e-5;             % thermal expansion, 1/K
MRHO(3,3)=1e-11;            % compressibility, 1/Pa
MRHO(3,4)=2400;             % melt density
MFLOW(3,1)=1;               % 1=power law (wet quartzite: Ranalli, 1995)
MFLOW(3,2)=3.2e-4;          % AD, 1/s/MPa^n
MFLOW(3,3)=2.3;             % n
MFLOW(3,4)=154;             % Ea, kJ/mol
MFLOW(3,5)=0;               % Va, cm^3
MMU(3,1)=1e+10;               % shear modulus, Pa
MPL(3,1)=1e+6;              % C0, Pa
MPL(3,2)=1e+6;              % C1, Pa
MPL(3,3)=0.15;              % sin(FI0) strong
MPL(3,4)=0.075;             % sin(FI1) weak
MPL(3,5)=0;                 % GAM0
MPL(3,6)=1;                 % GAM1
MCP(3,1)=1000;                % Cp, J/kg
MKT(3,1)=0.64;              % k0, W/m/K
MKT(3,2)=807;               % a, W/m
MHR(3,1)=2.0e-6;            % radiogenic heat production, W/m^3
MHE(3,1)=1;                 % shear and adiobatic heating efficiency

% -------------------------------------------------------------------------
% 4 = Upper oceanic crust (basalts)
% -------------------------------------------------------------------------
MRHO(4,1)=3000;             % standard density, kg/m^3
MRHO(4,2)=3e-5;             % thermal expansion, 1/K
MRHO(4,3)=1e-11;            % compressibility, 1/Pa
MRHO(4,4)=2900;             % melt density, kg/m^3
MFLOW(4,1)=1;               % 1=power law (wet quartzite: Ranalli, 1995)
MFLOW(4,2)=3.2e-4;          % AD, 1/s/MPa^n
MFLOW(4,3)=2.3;             % n
MFLOW(4,4)=154;             % Ea, kJ/mol
MFLOW(4,5)=0;               % Va, cm^3
MMU(4,1)=2.5e+10;             % shear modulus, Pa
MPL(4,1)=1e+6;              % C0, Pa
MPL(4,2)=1e+6;              % C1, Pa
MPL(4,3)=0.15;                 % sin(FI0), strong
MPL(4,4)=0.075;                 % sin(FI1), weak                
MPL(4,5)=0;                 % GAM0, strong  -- brittle strain accumulation
MPL(4,6)=1;                 % GAM1, weak    -- brittle strain accumulation
MCP(4,1)=1000;                % Cp, J/kg
MKT(4,1)=1.18;              % k0, W/m/K
MKT(4,2)=474;               % a, W/m
MHR(4,1)=2.5e-7;              % radiogenic heat production, W/m^3
MHE(4,1)=1;                 % shear and adiobatic heating efficiency

% -------------------------------------------------------------------------
% 5 = Lower oceanic crust (gabbro)
% -------------------------------------------------------------------------
MRHO(5,1)=3000;             % standard density, kg/m^3
MRHO(5,2)=3e-5;             % thermal expansion, 1/K
MRHO(5,3)=1e-11;            % compressibility, 1/Pa
MRHO(4,4)=2900;             % melt density, kg/m^3
MFLOW(5,1)=1;               % 1=power law (plagioclase An75: Ranalli, 1995)
MFLOW(5,2)=3.3e-4;          % AD, 1/s/MPa^n
MFLOW(5,3)=3.2;             % n
MFLOW(5,4)=238;             % Ea, kJ/mol
MFLOW(5,5)=0;               % Va, cm^3
MMU(5,1)=2.5e+10;             % shear modulus, Pa
MPL(5,1)=1e+6;              % C0, Pa
MPL(5,2)=1e+6;              % C1, Pa
MPL(5,3)=0.6;               % sin(FI0), strong
MPL(5,4)=0.3;               % sin(FI1), weak
MPL(5,5)=0;                 % GAM0
MPL(5,6)=1;                 % GAM1
MCP(5,1)=1000;                % Cp, J/kg
MKT(5,1)=1.18;              % k0, W/m/K
MKT(5,2)=474;               % a, W/m
MHR(5,1)=2.5e-7;              % radiogenic heat production, W/m^3
MHE(5,1)=1;                 % shear and adiobatic heating efficiency

% -------------------------------------------------------------------------
% 6 = Upper continental crust (granodiorite)
% -------------------------------------------------------------------------
MRHO(6,1)=2700;             % standard density, kg/m^3
MRHO(6,2)=3e-5;             % thermal expansion, 1/K
MRHO(6,3)=1e-11;            % compressibility, 1/Pa
MRHO(4,4)=2400;             % melt density, kg/m^3
MFLOW(6,1)=1;               % 1=power law (wet quartzite: Ranalli, 1995)
MFLOW(6,2)=3.2e-4;          % AD, 1/s/MPa^n
MFLOW(6,3)=2.3;             % n
MFLOW(6,4)=154;             % Ea, kJ/mol
MFLOW(6,5)=0;               % Va, cm^3
MMU(6,1)=1e+10;               % shear modulus, Pa
MPL(6,1)=1e+6;              % C0, Pa
MPL(6,2)=1e+6;              % C1, Pa
MPL(6,3)=0.15;                 % sin(FI0)
MPL(6,4)=0.075;                 % sin(FI1)
MPL(6,5)=0;                 % GAM0
MPL(6,6)=1;                 % GAM1
MCP(6,1)=1000;                % Cp, J/kg
MKT(6,1)=0.64;              % k0, W/m/K
MKT(6,2)=807;               % a, W/m
MHR(6,1)=1.0e-6;              % radiogenic heat production, W/m^3
MHE(6,1)=1;                 % shear and adiobatic heating efficiency

% -------------------------------------------------------------------------
% 7 = Lower continental crust (diorite)
% -------------------------------------------------------------------------
MRHO(7,1)=2900;             % standard density, kg/m^3
MRHO(7,2)=3e-5;             % thermal expansion, 1/K
MRHO(7,3)=1e-11;            % compressibility, 1/Pa
MRHO(4,4)=2400;             % melt density, kg/m^3
MFLOW(7,1)=1;               % 1=power law (plagioclase An75: Ranalli, 1995)
MFLOW(7,2)=3.3e-4;          % AD, 1/s/MPa^n
MFLOW(7,3)=3.2;             % n
MFLOW(7,4)=238;             % Ea, kJ/mol
MFLOW(7,5)=0;               % Va, cm^3
MMU(7,1)=2.5e+10;             % shear modulus, Pa
MPL(7,1)=1e+6;              % C0, Pa
MPL(7,2)=1e+6;              % C1, Pa
MPL(7,3)=0.15;               % sin(FI0)
MPL(7,4)=0.075;               % sin(FI1)
MPL(7,5)=0;                 % GAM0
MPL(7,6)=1;                 % GAM1
MCP(7,1)=1000;                % Cp, J/kg
MKT(7,1)=1.18;              % k0, W/m/K
MKT(7,2)=474;               % a, W/m
MHR(7,1)=5.0e-7;              % radiogenic heat production, W/m^3
MHE(7,1)=1;                 % shear and adiobatic heating efficiency

% -------------------------------------------------------------------------
% 8 = Oceanic Lithospheric mantle
% -------------------------------------------------------------------------
MRHO(8,1)=3300;             % standard density, kg/m^3
MRHO(8,2)=3e-5;             % thermal expansion, 1/K
MRHO(8,3)=1e-11;            % compressibility, 1/Pa
MFLOW(8,1)=1;               % 1=power law (dry olivine: Ranalli, 1995)
MFLOW(8,2)=2.5e+4;          % AD, 1/s/MPa^n
MFLOW(8,3)=3.5;             % n
MFLOW(8,4)=532;             % Ea, kJ/mol
MFLOW(8,5)=10;               % Va, cm^3
MMU(8,1)=6.7e+10;             % shear modulus, Pa
MPL(8,1)=1e+6;              % C0, Pa
MPL(8,2)=1e+6;              % C1, Pa
MPL(8,3)=0.6;               % sin(FI0)
MPL(8,4)=0.3;               % sin(FI1)
MPL(8,5)=0;                 % GAM0
MPL(8,6)=1;                 % GAM1
MCP(8,1)=1000;                % Cp, J/kg
MKT(8,1)=0.73;              % k0, W/m/K
MKT(8,2)=1293;              % a, W/m
MHR(8,1)=2.2e-8;              % radiogenic heat production, W/m^3
MHE(8,1)=1;                 % shear and adiobatic heating efficiency

% -------------------------------------------------------------------------
% 9 = Continetnal Lithospheric mantle -- dry mantle
% -------------------------------------------------------------------------
MRHO(9,1) = 3300;            % standard density, kg/m^3
MRHO(9,2) = 3e-5;            % thermal expansion, 1/K
MRHO(9,3) = 1e-11;           % compressibility, 1/Pa
MRHO(9,4) = 2900;            % melt density, kg/m^3   (Gerya and Mellick, 2010)
MRHO(9,5) = -15;             % depletion (compositional density constrast) kg/m^3
MFLOW(9,1)= 1;               % 1=power law (dry olivine: Ranalli, 1995)
MFLOW(9,2)= 2.5e+4;          % AD, 1/s/MPa^n
MFLOW(9,3)=3.5;              % n
MFLOW(9,4)=532;              % Ea, kJ/mol
MFLOW(9,5)=10;               % Va, cm^3
MMU(9)=6.7e+10;              % shear modulus, Pa
MPL(9,1)=10e+6;              % C0, Pa
MPL(9,2)=10e+6;              % C1, Pa
MPL(9,3)=0.6;                % sin(FI0)
MPL(9,4)=0.6;                % sin(FI1)
MPL(9,5)=0;                  % GAM0
MPL(9,6)=1;                  % GAM1
MCP(9,1)=1000;               % Cp, J/kg
MKT(9,1)=0.73;               % k0, W/m/K  (Gerya and Mellick, 2010)
MKT(9,2)=1293;               % a, W/m
MHR(9,1)=2.2e-8;             % radiogenic heat production, W/m^3
MHE(9,1)=1;                  % shear and adiobatic heating efficiency

% -------------------------------------------------------------------------
% 10 = Asthenospheric mantle -- dry
% -------------------------------------------------------------------------
MRHO(10,1)=3300;             % standard density, kg/m^3
MRHO(10,2)= 3e-5;             % thermal expansion, 1/K
MRHO(10,3)= 1e-11;            % compressibility, 1/Pa
MRHO(10,4)= 2900;             % melt density, kg/m^3
MFLOW(10,1)= 1;               % 1=power law (dry olivine: Ranalli, 1995)
MFLOW(10,2)= 2.5e+4;          % AD, 1/s/MPa^n
MFLOW(10,3)=3.5;             % n
MFLOW(10,4)=532;             % Ea, kJ/mol
MFLOW(10,5)=10;               % Va, cm^3
MMU(10,1)=6.7e+10;             % shear modulus, Pa
MPL(10,1)=10e+6;              % C0, Pa
MPL(10,2)=10e+6;              % C1, Pa
MPL(10,3)=0.6;               % sin(FI0)
MPL(10,4)=0.6;               % sin(FI1)
MPL(10,5)=0;                 % GAM0
MPL(10,6)=1;                 % GAM1
MCP(10,1)=1000;                % Cp, J/kg
MKT(10,1)=0.73;              % k0, W/m/K
MKT(10,2)=1293;              % a, W/m
MHR(10,1)=2.2e-8;              % radiogenic heat production, W/m^3
MHE(10,1)=1;                 % shear and adiobatic heating efficiency

% -------------------------------------------------------------------------
% 11 = Hydrated mantle in the intra-plate fracture zone (week zone)
% -------------------------------------------------------------------------
MRHO(11,1)=3300;             % standard density, kg/m^3
MRHO(11,2)=3e-5;             % thermal expansion, 1/K
MRHO(11,3)=1e-11;            % compressibility, 1/Pa
MRHO(11,4)=2900;             % melt density, kg/m^3   (Gerya and Mellick, 2010)
MRHO(11,5) = -15;            % depletion (compositional density constrast) kg/m^3
MFLOW(11,1)=1;               % 1=power law (wet olivine: Ranalli, 1995)
MFLOW(11,2)=2.0e+3;          % AD, 1/s/MPa^n
MFLOW(11,3)=4.0;             % n
MFLOW(11,4)=471;             % Ea, kJ/mol
MFLOW(11,5)=0;               % Va, cm^3
MMU(11,1)=6.7e+10;           % shear modulus, Pa
MPL(11,1)=1e+6;              % C0, Pa
MPL(11,2)=1e+6;              % C1, Pa
MPL(11,3)=0.1;               % sin(FI0)
MPL(11,4)=0.05;              % sin(FI1)
MPL(11,5)=0;                 % GAM0
MPL(11,6)=1;                 % GAM1
MCP(11,1)=1000;              % Cp, J/kg
MKT(11,1)=0.73;              % k0, W/m/K
MKT(11,2)=1293;              % a, W/m
MHR(11,1)=2.2e-8;            % radiogenic heat production, W/m^3
MHE(11,1)=1;                 % shear and adiobatic heating efficiency

% -------------------------------------------------------------------------
% 12 = Upper oceanic crust plateau (basalts)
% -------------------------------------------------------------------------
MRHO(12,1)=3000-150;             % standard density, kg/m^3
MRHO(12,2)=3e-5;             % thermal expansion, 1/K
MRHO(12,3)=1e-11;            % compressibility, 1/Pa
MFLOW(12,1)=1;               % 1=power law (wet quartzite: Ranalli, 1995)
MFLOW(12,2)=3.2e-4;          % AD, 1/s/MPa^n
MFLOW(12,3)=2.3;             % n
MFLOW(12,4)=154;             % Ea, kJ/mol
MFLOW(12,5)=0;               % Va, cm^3
MMU(12,1)=2.5e+10;             % shear modulus, Pa
MPL(12,1)=1e+6;              % C0, Pa
MPL(12,2)=1e+6;              % C1, Pa
MPL(12,3)=0.03;                 % sin(FI0)
MPL(12,4)=0.03;                 % sin(FI1)
MPL(12,5)=0;                 % GAM0
MPL(12,6)=1;                 % GAM1
MCP(12,1)=1000;                % Cp, J/kg
MKT(12,1)=1.18;              % k0, W/m/K
MKT(12,2)=474;               % a, W/m
MHR(12,1)=2.5e-7;              % radiogenic heat production, W/m^3
MHE(12,1)=1;                 % shear and adiobatic heating efficiency

% -------------------------------------------------------------------------
% 13 = Lower oceanic crust plateau
% -------------------------------------------------------------------------
MRHO(13,1)=3000-150;             % standard density, kg/m^3
MRHO(13,2)=3e-5;             % thermal expansion, 1/K
MRHO(13,3)=1e-11;            % compressibility, 1/Pa
MFLOW(13,1)=1;               % 1=power law (plagioclase An75: Ranalli, 1995)
MFLOW(13,2)=3.3e-4;          % AD, 1/s/MPa^n
MFLOW(13,3)=3.2;             % n
MFLOW(13,4)=238;             % Ea, kJ/mol
MFLOW(13,5)=0;               % Va, cm^3
MMU(13,1)=2.5e+10;             % shear modulus, Pa
MPL(13,1)=1e+6;              % C0, Pa
MPL(13,2)=1e+6;              % C1, Pa
MPL(13,3)=0.2;               % sin(FI0)
MPL(13,4)=0.2;               % sin(FI1)
MPL(13,5)=0;                 % GAM0
MPL(13,6)=1;                 % GAM1
MCP(13,1)=1000;                % Cp, J/kg
MKT(13,1)=1.18;              % k0, W/m/K
MKT(13,2)=474;               % a, W/m
MHR(13,1)=2.5e-7;              % radiogenic heat production, W/m^3
MHE(13,1)=1;                 % shear and adiobatic heating efficiency

% -------------------------------------------------------------------------
% 14 Plume material tracking equivlent to Asthenosphere = 6
% -------------------------------------------------------------------------
MRHO(14,1)=3300;             % standard density, kg/m^3
MRHO(14,2)=3e-5;             % thermal expansion, 1/K
MRHO(14,3)=1e-11;            % compressibility, 1/Pa
MRHO(14,4)=2900;             % melt density, kg/m^3
MRHO(14,5)= 0;             % compositional density constrast, kg/m^3
MFLOW(14,1)=1;               % 1=power law (dry olivine: Ranalli, 1995)
MFLOW(14,2)=2.5e+4;          % AD, 1/s/MPa^n
MFLOW(14,3)=3.5;             % n
MFLOW(14,4)=532;             % Ea, kJ/mol
MFLOW(14,5)=10;               % Va, cm^3
MMU(14,1)=6.7e+10;             % shear modulus, Pa
MPL(14,1)=10e+6;              % C0, Pa
MPL(14,2)=10e+6;              % C1, Pa
MPL(14,3)=0.6;               % sin(FI0)
MPL(14,4)=0.6;               % sin(FI1)
MPL(14,5)=0;                 % GAM0  : Accumulated strain regime
MPL(14,6)=1;                 % GAM1  : Accumulated strain regime
MCP(14,1)=1000;                % Cp, J/kg
MKT(14,1)=0.73;              % k0, W/m/K
MKT(14,2)=1293;              % a, W/m
MHR(14,1)=2.2e-8;              % radiogenic heat production, W/m^3% 10 Plume material tracking equivlent to Asthenosphere = 6
MHE(14,1)=1;                 % shear and adiobatic heating efficiency


% -------------------------------------------------------------------------
%% Defining lithological model constants
% -------------------------------------------------------------------------

% Sticky Layer thickness water or air (m)
sticky_layer = 20000;

% Water Level Tracking (m) set to empty if no water
% Water depth = sticky_layer - water_lev
water_lev = [];

% Sediment Thickness (m)
sed_thick = 2000;

% Continent level (m), use this to set continent above or at sea level
cont_lev = sticky_layer;

% Upper Continental Crust thickness (m)
upper_cont_thick = 23000 + sed_thick;  % Thickness includes sediment layer

% Lower Continental Crust thickness (m)
lower_cont_thick = 15000;

% Mantle lithosphere thickness (m)
mantle_lith_thick = 60000;

% Lithosphere thickness, including sediments (m)
litho_thick = mantle_lith_thick + upper_cont_thick + lower_cont_thick;

% Bottom of sediments
bot_sed = sticky_layer + sed_thick;

% Bottom of upper crust
bot_icrust1 = cont_lev + upper_cont_thick;

% Bottom of lower crust
bot_icrust2 = bot_icrust1 + lower_cont_thick;

% Continental lithosphere asthenosphere boundary (m) 
contLAB = sticky_layer + litho_thick;

% Plume head size and position (m)
rplume = 50000;
xplume = xsize/2;
yplume = 350000;

% Size and position of half-disk weak zone (m)
rweak = 5000;
xweak = xsize/2;
yweak = 60000; % brittle region of mantle lithosphere (model specific)

%--------------------------------------------------------------------------
%% Initial Temperature parameters
%--------------------------------------------------------------------------

% Weak zone temperature perturbation [k];
Tweak = 30;

% Plume Excess Temperature [K]
Texcess = 100;

% Thermal diffusivity of the mantle, m^2/s
kappa = 1e-6;    

% Convert deg-C to K
ttop = ttop - tabsolute;

% Temperature at bottom of model space
mantle_thick = ysize-sticky_layer;
tbottom = tpotential - tabsolute + tgrad*(mantle_thick);

% Temperature at bottom of continental lithosphere
tcontLAB = tpotential - tabsolute + tgrad*(contLAB - cont_lev);

%--------------------------------------------------------------------------
%% Generate initial strength envelopes
%--------------------------------------------------------------------------

% Temperature Profile
SE.z = linspace(0,mantle_thick,(1+mantle_thick/1000))';

% Mantle temperature (adiabatic)
SE.T = tbottom - tgrad*(mantle_thick-SE.z);
k = SE.z <= litho_thick;
temp = SE.T(k);
z = SE.z(k);

% Linear temperature from ttop to tcontLAB
t_slope = (temp(end) - ttop) / litho_thick;
t_intercept = temp(end) - z(end)*t_slope;
SE.T(k) = t_slope*SE.z(k) + t_intercept;

% Continental Lithosphere
SE = strengthprofile(sed_thick,MRHO(3,:),MFLOW(3,:),MPL(3,:),...
    upper_cont_thick-sed_thick,MRHO(6,:),MFLOW(6,:),MPL(6,:),...
    lower_cont_thick,MRHO(7,:),MFLOW(7,:),MPL(7,:),...
    mantle_thick,MRHO(10,:),MFLOW(10,:),MPL(10,:),SE);

save([outpath,'/litho_strength.mat'],'SE'); 
% -------------------------------------------------------------------------
%% Defining lithological structure of the model
% -------------------------------------------------------------------------

% Marker counter
mm1=0;
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
    if(MY(mm1)>sticky_layer && MY(mm1)<= bot_sed)
        MI(mm1)=3;
    end
    
    % Upper continental crust
    if(MY(mm1)>bot_sed && MY(mm1)<=bot_icrust1)
        MI(mm1)=6;
    end
    
    % Lower continental crust
    if(MY(mm1)>bot_icrust1 && MY(mm1)<=bot_icrust2)
        MI(mm1)=7;
    end
    
    % Continetal lab
    if(MY(mm1)>bot_icrust2 && MY(mm1)<=contLAB)
        MI(mm1) = 9;
    end
    
    %--------------------------------------------------------
    % Plume Head
    %--------------------------------------------------------
      
    dx=MX(mm1) - xplume;
    dy=MY(mm1) - yplume;
    dr_plume = (dx^2+dy^2)^0.5;
    if(dr_plume <= rplume)
        MI(mm1)=14;
    end
    
    
    %--------------------------------------------------------
    % Initial temperature structure
    %--------------------------------------------------------
    
    % Adiabatic Temperature Gradient in the Mantle
    MTK(mm1)=tbottom-tgrad*(ysize-MY(mm1));
    
    % Sticky air or sticky water
    if(MI(mm1) == 1 || MI(mm1) == 2)
        MTK(mm1) = ttop;
    end
    
    % Continental geotherm is linear from surface to LAB
    if( MY(mm1) > cont_lev && MY(mm1) <= contLAB)
        MTK(mm1) = ttop + (tcontLAB-ttop)*(MY(mm1)-cont_lev)/(contLAB-cont_lev);
    end
    
    %--------------------------------------------------------
    % Thermal Weak zone
    %--------------------------------------------------------
    
    % Thermal perturbation in the base of the lithosphere
    % (to start extension in the middle of the model)
    dx = MX(mm1) - xweak;
    dy = MY(mm1) - yweak;
    dr = (dx^2+dy^2)^0.5;
    if (MY(mm1) <= yweak && dr <= rweak)
        MI(mm1) = 11;
%         MTK(mm1) = MTK(mm1) + Tweak;%*(1-dr/rweak);
    end
    
    %--------------------------------------------------------
    % Plume Temperature
    %--------------------------------------------------------
    
    if (MI(mm1) == 14)
        MTK(mm1) = MTK(mm1) + Texcess*(1-dr_plume/rplume);  % Gaussian temperature
    end
    
end

% Save Number of markers marknum=mm1
marknum = mm1;

% -------------------------------------------------------------------------
%% Boundary Conditions -- Velocity
% -------------------------------------------------------------------------
btop = zeros(xnum+1,4);
bbottom = zeros(xnum+1,4);
bleft = zeros(ynum+1,4);
bright = zeros(ynum+1,4);

% Upper, Lower boundaries: Free slip + Prescribed inward velocity (vertical shortening)
% Velocities are determined according to equations by Laoi and Gerya,
% Tectonophysics, 2014

% Upper boundary: Free slip
% vx(1,j)=btop(j,1)+vx(2,j)*btop(j,2)
btop(:,1)=0;
btop(:,2)=1;
% vy(1,j)=btop(j,3)+vy(2,j)*btop(j,4)
btop(:,3)=0;
btop(:,4)= sticky_layer*(vxright-vxleft)/xsize;

% Lower boundary: Free slip
% vx(ynum+1,j)=bbottom(j,1)+vx(ynum,j)*bbottom(j,2)
bbottom(:,1)= 0;
bbottom(:,2)= 1;
% vy(ynum,j)=bbottom(j,3)+vy(ynum-1,j)*bbottom(j,4)
bbottom(:,3)= (vxright-vxleft)*(sticky_layer-ysize)/xsize;
bbottom(:,4)= 0;

% % Lower boundary: Vx Free slip, Vy(zext=1000 km) = 0 [Permeable]
% dz = gridy(ynum)-gridy(ynum-1);
% zext = 1000e3;
% % vx(ynum+1,j)=bbottom(j,1)+vx(ynum,j)*bbottom(j,2)
% bbottom(:,1)=0;
% bbottom(:,2)=1;
% % vy(ynum,j)=bbottom(j,3)+vy(ynum-1,j)*bbottom(j,4)
% bbottom(:,3)=0;
% bbottom(:,4)=zext/(zext+dz);

% Left boundary: Free slip, Vy = 25 mm/yr extension
% vx(i,1)=bleft(i,1)+vx(i,2)*bleft(i,2)
bleft(:,1)=vxleft;
bleft(:,2)=0;
% vy(i,1)=bleft(i,3)+vy(i,2)*bleft(i,4)
bleft(:,3)=0;
bleft(:,4)=1;

% Right boundary: Free slip, Vy = 25 mm/yr extension
% vx(i,xnum)=bright(i,1)+vx(i,xnum-1)*bbright(i,2)
bright(:,1)=vxright;
bright(:,2)=0;
% vy(i,xnum+1)=bright(i,3)+vx(i,xnum)*bbright(i,4)
bright(:,3)=0;
bright(:,4)=1;

% -------------------------------------------------------------------------
%% Boundary Conditions -- Thermal
% -------------------------------------------------------------------------
temp_condition = false;

btopt = zeros(xnum,2);
bbottomt = btopt;
bleftt = zeros(ynum,2);
brightt = bleftt;

%---------------------------------------------
% Upper Boundary fixed
%---------------------------------------------
% tk(1,j)=btopt(j,1)+tk(2,j)*btop(j,2)
btopt(:,1) = ttop;
btopt(:,2) = 0;

%---------------------------------------------
% Lower Boundary fixed
%---------------------------------------------
bbottomt(:,1) = tbottom;
bbottomt(:,2) = 0;


% Add temperature anomaly corresponding to plume
temp_condition = true;
i = gridx >= xplume-rplume & gridx <= xplume+rplume;
bbottomt(i,1) = bbottomt(i,1) + Texcess;

% % Lower Boundary external fixed (grids can't change)
% % tk(ynum,j)=bbottomt(j,1)+tk(ynum-1,j)*bbottomt(j,2)
% % bbottomt(:,1) = tbottom;
% dz = gridy(ynum)-gridy(ynum-1);
% zext = 1000e3;
% Text = tpotential - tabsolute + tgrad*(zext+gridy(ynum)-sticky_layer);
% % 
% bbottomt(:,1) = Text*dz/(zext+dz);
% bbottomt(:,2) = zext/(zext+dz);

%---------------------------------------------
% Left, Right boundaries: symmetry (zero flux)
%---------------------------------------------
for i=1:1:ynum
    % Left boundary
    % tk(i,1)=bleftt(i,1)+bleftt(i,2)*tk(i,2);
    bleftt(i,1)=0;
    bleftt(i,2)=1;
    % Right boundary
    % tk(i,xnum)=brightt(i,1)+brightt(i,2)*tk(i,xnum-1);
    brightt(i,1)=0;
    brightt(i,2)=1;    
end

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