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
MatProp{15} = 'Eclogite';

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
% MKT = thermal conductivity (W/m/K): k=k0+a/(T+77)*(1+b*P)
% MHR = radiogenic heat production (W/m^3)
% MHE = shear and adiobatic heating efficiency:  take on values 0 or 1
% MPH = phase change properties, clapyron slope, delta_rho

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
MPL(1,7)=0;                 % Maximum yield strengh
MCP(1,1)=3000;              %Cp, J/kg
MKT(1,1)=300;               % k0, W/m/K
MKT(1,2)=0;                 % a, W/m
MKT(1,3)=0;                 % Pressure dependance of thermal conductivity b, W/m/Pa
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
MPL(2,7)=0;                 % Maximum yield strengh
MCP(2,1)=3000;                % Cp, J/kg
MKT(2,1)=300;               % k0, W/m/K
MKT(2,2)=0;                 % a, W/m
MKT(2,3)=0;                 % Pressure dependance of thermal conductivity b, W/m/Pa
MHR(2,1)=0;                 % radiogenic heat production, W/m^3
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
MFLOW(3,5)=0;               % Va, cm^3 ===> *1e-6 for J/Pa, or *0.1 for J/bar
MMU(3,1)=1e+10;               % shear modulus, Pa
MPL(3,1)=10e+6;              % C0, Pa
MPL(3,2)=3e+6;              % C1, Pa
MPL(3,3)=0;              % sin(FI0) strong
MPL(3,4)=0;             % sin(FI1) weak
MPL(3,5)=0;                 % GAM0
MPL(3,6)=1;                 % GAM1
MPL(3,7)=200e6;               % Maximum yield strength, Pa
MCP(3,1)=1000;                % Cp, J/kg
MKT(3,1)=0.64;              % k0, W/m/K
MKT(3,2)=807;               % a, W/m
MKT(3,3)=0;                 % Pressure dependance of thermal conductivity b, W/m/Pa
MHR(3,1)=1.0e-6;            % radiogenic heat production, W/m^3
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
MFLOW(4,5)=8;               % Va, cm^3 ===> *1e-6 for J/Pa, or *0.1 for J/bar
MMU(4,1)=2.5e+10;             % shear modulus, Pa
MPL(4,1)=1e+6;              % C0, Pa
MPL(4,2)=1e+6;              % C1, Pa
MPL(4,3)=0.1;                 % sin(FI0), strong
MPL(4,4)=0.1;                 % sin(FI1), weak                
MPL(4,5)=0;                 % GAM0, strong  -- brittle strain accumulation
MPL(4,6)=1;                 % GAM1, weak    -- brittle strain accumulation
MPL(4,7)=200e6;             % Maximum yield strength, Pa
MCP(4,1)=1000;                % Cp, J/kg
MKT(4,1)=1.18;              % k0, W/m/K
MKT(4,2)=474;               % a, W/m
MKT(4,3)=0;                 % Pressure dependance of thermal conductivity b, W/m/Pa
MHR(4,1)=2.5e-7;              % radiogenic heat production, W/m^3
MHE(4,1)=1;                 % shear and adiobatic heating efficiency

% -------------------------------------------------------------------------
% 5 = Lower oceanic crust (gabbro)
% -------------------------------------------------------------------------
MRHO(5,1)=3000;             % standard density, kg/m^3
MRHO(5,2)=3e-5;             % thermal expansion, 1/K
MRHO(5,3)=1e-11;            % compressibility, 1/Pa
MRHO(5,4)=2900;             % melt density, kg/m^3
MFLOW(5,1)=1;               % 1=power law (plagioclase An75: Ranalli, 1995)
MFLOW(5,2)=3.3e-4;          % AD, 1/s/MPa^n
MFLOW(5,3)=3.2;             % n
MFLOW(5,4)=238;             % Ea, kJ/mol
MFLOW(5,5)=8;               % Va, cm^3  ===> *1e-6 for J/Pa, or *0.1 for J/bar
MMU(5,1)=2.5e+10;             % shear modulus, Pa
MPL(5,1)=1e+6;              % C0, Pa
MPL(5,2)=1e+6;              % C1, Pa
MPL(5,3)=0.1;               % sin(FI0), strong
MPL(5,4)=0.1;               % sin(FI1), weak
MPL(5,5)=8;                 % GAM0
MPL(5,6)=1;                 % GAM1
MPL(5,7)=200e6;               % Maximum yield strength, Pa
MCP(5,1)=1000;                % Cp, J/kg
MKT(5,1)=1.18;              % k0, W/m/K
MKT(5,2)=474;               % a, W/m
MKT(5,3)=0;                 % Pressure dependance of thermal conductivity b, W/m/Pa
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
MFLOW(6,5)=0;               % Va, cm^3 ===> *1e-6 for J/Pa, or *0.1 for J/bar
MMU(6,1)=1e+10;               % shear modulus, Pa
MPL(6,1)=10e+6;              % C0, Pa
MPL(6,2)=3e+6;              % C1, Pa
MPL(6,3)=0.6;                 % sin(FI0)
MPL(6,4)=0.3;                 % sin(FI1)
MPL(6,5)=0;                 % GAM0
MPL(6,6)=0.25;                 % GAM1
MPL(6,7)=200e6;               % Maximum yield strength, Pa
MCP(6,1)=1000;                % Cp, J/kg
MKT(6,1)=0.64;              % k0, W/m/K
MKT(6,2)=807;               % a, W/m
MKT(6,3)=0;                 % Pressure dependance of thermal conductivity b, W/m/Pa
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
MFLOW(7,5)=0;               % Va, cm^3 ===> *1e-6 for J/Pa, or *0.1 for J/bar
MMU(7,1)=2.5e+10;             % shear modulus, Pa
MPL(7,1)=10e+6;              % C0, Pa
MPL(7,2)=3e+6;              % C1, Pa
MPL(7,3)=0.6;               % sin(FI0)
MPL(7,4)=0.3;               % sin(FI1)
MPL(7,5)=0;                 % GAM0
MPL(7,6)=0.25;                 % GAM1
MPL(7,7)=200e6;               % Maximum yield strength, Pa
MCP(7,1)=1000;                % Cp, J/kg
MKT(7,1)=1.18;              % k0, W/m/K
MKT(7,2)=474;               % a, W/m
MKT(7,3)=0;                 % Pressure dependance of thermal conductivity b, W/m/Pa
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
MFLOW(8,5)=12;               % Va, cm^3 ===> *1e-6 for J/Pa, or *0.1 for J/bar
MMU(8,1)=6.7e+10;             % shear modulus, Pa
MPL(8,1)=1e+6;              % C0, Pa
MPL(8,2)=1e+6;              % C1, Pa
MPL(8,3)=0.15;               % sin(FI0)
MPL(8,4)=0.15;               % sin(FI1)
MPL(8,5)=0;                 % GAM0
MPL(8,6)=1;                 % GAM1
MPL(8,7)=200e6;               % Maximum yield strength, Pa
MCP(8,1)=1000;                % Cp, J/kg
MKT(8,1)=1.73;              % k0, W/m/K
MKT(8,2)=1293;              % a, W/m
MKT(8,3)=4e-11;             % Pressure dependance of thermal conductivity b, W/m/Pa
MHR(8,1)=2.2e-8;              % radiogenic heat production, W/m^3
MHE(8,1)=1;                 % shear and adiobatic heating efficiency

% -------------------------------------------------------------------------
% 9 = Continental Lithospheric mantle -- dry mantle
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
MFLOW(9,5)=10;               % Va, cm^3 ===> *1e-6 for J/Pa, or *0.1 for J/bar
MMU(9)=6.7e+10;              % shear modulus, Pa
MPL(9,1)=10e+6;              % C0, Pa
MPL(9,2)=3e+6;              % C1, Pa
MPL(9,3)=0.6;                % sin(FI0)
MPL(9,4)=0.3;                % sin(FI1)
MPL(9,5)=0;                  % GAM0
MPL(9,6)=1;                  % GAM1
MPL(9,7)=200e6;               % Maximum yield strength, Pa
MCP(9,1)=1000;               % Cp, J/kg
MKT(9,1)=0.73;               % k0, W/m/K  (Gerya and Mellick, 2010)
MKT(9,2)=1293;               % a, W/m
MKT(9,3)=4e-11;             % Pressure dependance of thermal conductivity b, W/m/Pa
MHR(9,1)=4e-9;             % radiogenic heat production, W/m^3
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
MFLOW(10,5)=12;               % Va, cm^3 ===> *1e-6 for J/Pa, or *0.1 for J/bar
MMU(10,1)=6.7e+10;             % shear modulus, Pa
MPL(10,1)=3e+6;              % C0, Pa
MPL(10,2)=3e+6;              % C1, Pa
MPL(10,3)=0.1;               % sin(FI0)
MPL(10,4)=0.0;               % sin(FI1)
MPL(10,5)=0;                 % GAM0
MPL(10,6)=1;                 % GAM1
MPL(10,7)=200e6;               % Maximum yield strength, Pa
MCP(10,1)=1000;                % Cp, J/kg
MKT(10,1)=1.73;              % k0, W/m/K
MKT(10,2)=1293;              % a, W/m
MKT(10,3)=4e-11;             % Pressure dependance of thermal conductivity b, W/m/Pa
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
MFLOW(11,5)=5;               % Va, cm^3 ===> *1e-6 for J/Pa, or *0.1 for J/bar
MMU(11,1)=6.7e+10;           % shear modulus, Pa
MPL(11,1)=1e+6;              % C0, Pa
MPL(11,2)=1e+6;              % C1, Pa
MPL(11,3)=0;               % sin(FI0)
MPL(11,4)=0;              % sin(FI1)
MPL(11,5)=0;                 % GAM0
MPL(11,6)=1;                 % GAM1
MPL(11,7)=200e6;               % Maximum yield strength, Pa
MCP(11,1)=1000;              % Cp, J/kg
MKT(11,1)=0.73;              % k0, W/m/K
MKT(11,2)=1293;              % a, W/m
MKT(11,3)=4e-11;             % Pressure dependance of thermal conductivity b, W/m/Pa
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
MFLOW(12,5)=0;               % Va, cm^3 ===> *1e-6 for J/Pa, or *0.1 for J/bar
MMU(12,1)=2.5e+10;             % shear modulus, Pa
MPL(12,1)=1e+6;              % C0, Pa
MPL(12,2)=1e+6;              % C1, Pa
MPL(12,3)=0.03;                 % sin(FI0)
MPL(12,4)=0.03;                 % sin(FI1)
MPL(12,5)=0;                 % GAM0
MPL(12,6)=1;                 % GAM1
MPL(12,7)=200e6;               % Maximum yield strength, Pa
MCP(12,1)=1000;                % Cp, J/kg
MKT(12,1)=1.18;              % k0, W/m/K
MKT(12,2)=474;               % a, W/m
MKT(12,3)=0;                 % Pressure dependance of thermal conductivity b, W/m/Pa
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
MFLOW(13,5)=0;               % Va, cm^3 ===> *1e-6 for J/Pa, or *0.1 for J/bar
MMU(13,1)=2.5e+10;             % shear modulus, Pa
MPL(13,1)=1e+6;              % C0, Pa
MPL(13,2)=1e+6;              % C1, Pa
MPL(13,3)=0.2;               % sin(FI0)
MPL(13,4)=0.2;               % sin(FI1)
MPL(13,5)=0;                 % GAM0
MPL(13,6)=1;                 % GAM1
MPL(13,7)=200e6;               % Maximum yield strength, Pa
MCP(13,1)=1000;                % Cp, J/kg
MKT(13,1)=1.18;              % k0, W/m/K
MKT(13,2)=474;               % a, W/m
MKT(13,3)=0;                 % Pressure dependance of thermal conductivity b, W/m/Pa
MHR(13,1)=2.5e-7;              % radiogenic heat production, W/m^3
MHE(13,1)=1;                 % shear and adiobatic heating efficiency

% -------------------------------------------------------------------------
% 14 Plume material tracking equivlent to Asthenosphere = 6
% -------------------------------------------------------------------------
MRHO(14,1)=3300;             % standard density, kg/m^3
MRHO(14,2)=3e-5;             % thermal expansion, 1/K
MRHO(14,3)=1e-11;            % compressibility, 1/Pa
MRHO(14,4)=2900;             % melt density, kg/m^3
MRHO(14,5)= -10;             % compositional density constrast, kg/m^3
MFLOW(14,1)=1;               % 1=power law (dry olivine: Ranalli, 1995)
MFLOW(14,2)=2.5e+4;          % AD, 1/s/MPa^n
MFLOW(14,3)=3.5;             % n
MFLOW(14,4)=532;             % Ea, kJ/mol
MFLOW(14,5)=10;               % Va, cm^3 ===> *1e-6 for J/Pa, or *0.1 for J/bar
MMU(14,1)=6.7e+10;             % shear modulus, Pa
MPL(14,1)=10e+6;              % C0, Pa
MPL(14,2)=10e+6;              % C1, Pa
MPL(14,3)=0.6;               % sin(FI0)
MPL(14,4)=0.6;               % sin(FI1)
MPL(14,5)=0;                 % GAM0  : Accumulated strain regime
MPL(14,6)=1;                 % GAM1  : Accumulated strain regime
MPL(14,7)=200e6;               % Maximum yield strength, Pa
MCP(14,1)=1000;                % Cp, J/kg
MKT(14,1)=0.73;              % k0, W/m/K
MKT(14,2)=1293;              % a, W/m
MKT(14,3)=4e-11;             % Pressure dependance of thermal conductivity b, W/m/Pa
MHR(14,1)=2.2e-8;              % radiogenic heat production, W/m^3% 10 Plume material tracking equivlent to Asthenosphere = 6
MHE(14,1)=1;                 % shear and adiobatic heating efficiency

% -------------------------------------------------------------------------
% 15 = Eclogite (Use Gabbro Rheology instead of basalt) -- not used
% -------------------------------------------------------------------------
MRHO(15,1)=3400;             % standard density, kg/m^3
MRHO(15,2)=3e-5;             % thermal expansion, 1/K
MRHO(15,3)=1e-11;            % compressibility, 1/Pa
MRHO(15,4)=2900;             % melt density, kg/m^3
MFLOW(15,1)=1;               % 1=power law (plagioclase An75: Ranalli, 1995)
MFLOW(15,2)=3.3e-4;          % AD, 1/s/MPa^n
MFLOW(15,3)=3.2;             % n
MFLOW(15,4)=480;             % Ea, kJ/mol
MFLOW(15,5)=0;               % Va, cm^3 ===> *1e-6 for J/Pa, or *0.1 for J/bar
MMU(15,1)=2.5e+10;           % shear modulus, Pa
MPL(15,1)=1e+6;              % C0, Pa
MPL(15,2)=1e+6;              % C1, Pa
MPL(15,3)=0.6;               % sin(FI0), strong
MPL(15,4)=0.3;               % sin(FI1), weak
MPL(15,5)=0;                 % GAM0
MPL(15,6)=1;                 % GAM1
MPL(15,7)=200e6;               % Maximum yield strength, Pa
MCP(15,1)=1000;              % Cp, J/kg
MKT(15,1)=1.18;              % k0, W/m/K
MKT(15,2)=474;               % a, W/m
MKT(15,3)=0;                 % Pressure dependance of thermal conductivity b, W/m/Pa
MHR(15,1)=2.5e-7;            % radiogenic heat production, W/m^3
MHE(15,1)=1;                 % shear and adiobatic heating efficiency

% -------------------------------------------------------------------------
%% Defining Phase Transformations
% (1) Basalt to eclogite
% (2) 410 km   see Faccenda et al (2017) for parameters
% (3) 660 km
% -------------------------------------------------------------------------
MPHASE = zeros(9,4);

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
MPHASE(5,1:2) = 100;              % Density Change (kg/m^3)

% 410 km
MPHASE(4,3) = 3;                  % Clapeyron slope (MPa/K)
MPHASE(5,3) = 206;                % Density change (kg/m^3)
MPHASE(6,3) = 13.6;               % Reference phase change pressure (GPa)
MPHASE(7,3) = 1810;               % Reference phase change temperature (K)
MPHASE(8,3) = 1;                  % Viscosity increase factor
MPHASE(9,3) = 0.1;                % Phase transition range for hyperbolic (0.1 GPa = 20 km)

% 660 km
MPHASE(4,4) = -2.5;               % Clapeyron slope (MPa/K)
MPHASE(5,4) = 322;                % Density change (kg/m^3)
MPHASE(6,4) = 23.3;              % Reference phase change pressure (GPa)
MPHASE(7,4) = 1940;               % Reference phase change temperature (K)
MPHASE(8,4) = 30;                 % Viscosity increase factor
MPHASE(9,4) = 0.1;                % Phase transition range for hyperbolic (0.1 GPa = 20 km)

% -------------------------------------------------------------------------
%% Defining lithological temperature model constants
% -------------------------------------------------------------------------

% Sticky Layer thickness water or air (m)
sticky_layer = 20000;

% Water Level Tracking (m) set to empty if no water
% Water depth = sticky_layer - water_lev
water_lev = 14000;

% For plume perturbation
add_plume = false;

%--------------------------------------------------------------------------
%% Initial Temperature parameters
%--------------------------------------------------------------------------

% Convert deg-C to K
ttop = ttop - tabsolute;

% Temperature at bottom of model space
mantle_thick = ysize-sticky_layer;
tbottom = tpotential - tabsolute + tgrad*mantle_thick;

% Thermal lithosphere
ThermLAB = 1100 - tabsolute;

% Compute subduction temperatures
mm1 = MY > sticky_layer & (MY-sticky_layer) <= 200e3;
MTK(mm1) = OceanOceanSubThermal(MX(mm1)/1000,(MY(mm1)-sticky_layer)/1000,tpotential) - tabsolute;

% -------------------------------------------------------------------------
%% Defining lithological structure of the model
% -------------------------------------------------------------------------

% Weak zone polygon
theta = 1:65;
x1 = 1800e3 + 200e3*cosd(90-theta);
x2 = 1800e3 + 210e3*cosd(90-theta);
y1 = 200e3 - 200e3*sind(90-theta) + sticky_layer;
y2 = 200e3 - 210e3*sind(90-theta) + sticky_layer;
k = y1 < sticky_layer+7e3;
x1(k) = [];
y1(k) = [];
k = y2 < sticky_layer+7e3;
x2(k) = [];
y2(k) = [];

xp = [x1 fliplr(x2)];
yp = [y1 fliplr(y2)];
        
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
      
    %--------------------------------------------------------
    % Thermal models
    %--------------------------------------------------------
    
    % Adiabatic Temperature Gradient in the Mantle
    if (MY(mm1)-sticky_layer) > 120e3
        MTK(mm1)=tbottom-tgrad*(ysize-MY(mm1));
    end
    
    % Sticky air or sticky water
    if(MI(mm1) == 1 || MI(mm1) == 2)
        MTK(mm1) = ttop;
    end
    
    %--------------------------------------------------------
    % Oceanic Plates defined by temperature
    %--------------------------------------------------------
    
    % Mantle lithosphere
    if(MY(mm1)>sticky_layer && MTK(mm1) <= ThermLAB)
        MI(mm1) = 8;
    end    
    
    % Basalt Crust
    if MY(mm1)>sticky_layer && MY(mm1) <= sticky_layer + 9e3
        MI(mm1) = 4;
    end
    
    % Sediment
    if MY(mm1)>sticky_layer && MY(mm1) <= sticky_layer + 2e3
        MI(mm1) = 3;
    end
    
end
% Create weak zone
[in,on]=inpolygon(MX,MY,xp,yp);
MI(in|on) = 11;


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