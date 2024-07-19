% Mantle Convection

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
MPL(2,7)=0;                 % Maximum yield strengh
MCP(2,1)=3000;                % Cp, J/kg
MKT(2,1)=300;               % k0, W/m/K
MKT(2,2)=0;                 % a, W/m
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
MPL(3,7)=1000e6;               % Maximum yield strength, Pa
MCP(3,1)=1000;                % Cp, J/kg
MKT(3,1)=0.64;              % k0, W/m/K
MKT(3,2)=807;               % a, W/m
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
MPL(4,7)=1000e6;             % Maximum yield strength, Pa
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
MPL(5,7)=1000e6;               % Maximum yield strength, Pa
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
MPL(6,7)=1000e6;               % Maximum yield strength, Pa
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
MFLOW(7,5)=0;               % Va, cm^3 ===> *1e-6 for J/Pa, or *0.1 for J/bar
MMU(7,1)=2.5e+10;             % shear modulus, Pa
MPL(7,1)=10e+6;              % C0, Pa
MPL(7,2)=3e+6;              % C1, Pa
MPL(7,3)=0.6;               % sin(FI0)
MPL(7,4)=0.3;               % sin(FI1)
MPL(7,5)=0;                 % GAM0
MPL(7,6)=0.25;                 % GAM1
MPL(7,7)=1000e6;               % Maximum yield strength, Pa
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
MFLOW(8,5)=8;               % Va, cm^3 ===> *1e-6 for J/Pa, or *0.1 for J/bar
MMU(8,1)=6.7e+10;             % shear modulus, Pa
MPL(8,1)=1e+6;              % C0, Pa
MPL(8,2)=1e+6;              % C1, Pa
MPL(8,3)=0.15;               % sin(FI0)
MPL(8,4)=0.15;               % sin(FI1)
MPL(8,5)=0;                 % GAM0
MPL(8,6)=1;                 % GAM1
MPL(8,7)=1000e6;               % Maximum yield strength, Pa
MCP(8,1)=1000;                % Cp, J/kg
MKT(8,1)=0.73;              % k0, W/m/K
MKT(8,2)=1293;              % a, W/m
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
MPL(9,7)=1000e6;               % Maximum yield strength, Pa
MCP(9,1)=1000;               % Cp, J/kg
MKT(9,1)=0.73;               % k0, W/m/K  (Gerya and Mellick, 2010)
MKT(9,2)=1293;               % a, W/m
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
MFLOW(10,5)=5;               % Va, cm^3 ===> *1e-6 for J/Pa, or *0.1 for J/bar
MMU(10,1)=6.7e+10;             % shear modulus, Pa
MPL(10,1)=1e+6;              % C0, Pa
MPL(10,2)=1e+6;              % C1, Pa
MPL(10,3)=0.1;               % sin(FI0)
MPL(10,4)=0.0;               % sin(FI1)
MPL(10,5)=0;                 % GAM0
MPL(10,6)=1;                 % GAM1
MPL(10,7)=1000e6;               % Maximum yield strength, Pa
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
MFLOW(11,5)=0;               % Va, cm^3 ===> *1e-6 for J/Pa, or *0.1 for J/bar
MMU(11,1)=6.7e+10;           % shear modulus, Pa
MPL(11,1)=1e+6;              % C0, Pa
MPL(11,2)=1e+6;              % C1, Pa
MPL(11,3)=0;               % sin(FI0)
MPL(11,4)=0;              % sin(FI1)
MPL(11,5)=0;                 % GAM0
MPL(11,6)=1;                 % GAM1
MPL(11,7)=1000e6;               % Maximum yield strength, Pa
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
MFLOW(12,5)=0;               % Va, cm^3 ===> *1e-6 for J/Pa, or *0.1 for J/bar
MMU(12,1)=2.5e+10;             % shear modulus, Pa
MPL(12,1)=1e+6;              % C0, Pa
MPL(12,2)=1e+6;              % C1, Pa
MPL(12,3)=0.03;                 % sin(FI0)
MPL(12,4)=0.03;                 % sin(FI1)
MPL(12,5)=0;                 % GAM0
MPL(12,6)=1;                 % GAM1
MPL(12,7)=1000e6;               % Maximum yield strength, Pa
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
MFLOW(13,5)=0;               % Va, cm^3 ===> *1e-6 for J/Pa, or *0.1 for J/bar
MMU(13,1)=2.5e+10;             % shear modulus, Pa
MPL(13,1)=1e+6;              % C0, Pa
MPL(13,2)=1e+6;              % C1, Pa
MPL(13,3)=0.2;               % sin(FI0)
MPL(13,4)=0.2;               % sin(FI1)
MPL(13,5)=0;                 % GAM0
MPL(13,6)=1;                 % GAM1
MPL(13,7)=1000e6;               % Maximum yield strength, Pa
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
MPL(14,7)=1000e6;               % Maximum yield strength, Pa
MCP(14,1)=1000;                % Cp, J/kg
MKT(14,1)=0.73;              % k0, W/m/K
MKT(14,2)=1293;              % a, W/m
MHR(14,1)=2.2e-8;              % radiogenic heat production, W/m^3% 10 Plume material tracking equivlent to Asthenosphere = 6
MHE(14,1)=1;                 % shear and adiobatic heating efficiency

% -------------------------------------------------------------------------
% 15 = Eclogite (Use Gabbro Rheology instead of basalt)
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
MPL(15,7)=1000e6;               % Maximum yield strength, Pa
MCP(15,1)=1000;              % Cp, J/kg
MKT(15,1)=1.18;              % k0, W/m/K
MKT(15,2)=474;               % a, W/m
MHR(15,1)=2.5e-7;            % radiogenic heat production, W/m^3
MHE(15,1)=1;                 % shear and adiobatic heating efficiency

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
MPHASE(5,3) = 206;              % Density change (kg/m^3)
MPHASE(6,3) = 13.5;             % Reference phase change pressure (GPa)
MPHASE(7,3) = 1810;             % Reference phase change temperature (K)
MPHASE(8,3) = 1;                % Viscosity increase factor

% 660 km
MPHASE(4,4) = -1.3;             % Clapeyron slope (MPa/K)
MPHASE(5,4) = 322;              % Density change (kg/m^3)
MPHASE(6,4) = 23;               % Reference phase change pressure (GPa)
MPHASE(7,4) = 1940;             % Reference phase change temperature (K)
MPHASE(8,4) = 1;                % Viscosity increase factor

% -------------------------------------------------------------------------
%% Defining lithological temperature model constants
% -------------------------------------------------------------------------

% Sticky Layer thickness water or air (m)
sticky_layer = 0;

% Water Level Tracking (m) set to empty if no water
% Water depth = sticky_layer - water_lev
water_lev = [];

% Sediment Thickness (m)
sed_thick = 0;

% Thermal diffusivity of the mantle, m^2/s
kappa = 1e-6;    

% Oceanic lithosphere thickness (1200oC isotherm) (m)
isoLAB = 1100-tabsolute;

% Age of Oceanic Plate (s)
plate_age_old = 60e6*yr2sec;  % Subducting

% Age of Oceanic Plate (s)
plate_age_young = 30e6*yr2sec; % Overriding

%--------------------------------------------------------------------------
%% Initial Temperature parameters
%--------------------------------------------------------------------------

% Convert deg-C to K
ttop = ttop - tabsolute;
tbottom = 3500;

% Temperature at bottom of model space
mantle_thick = ysize-sticky_layer;

% Upper depth
udepth = 200e3;

% Half-space cooling geotherm
z = linspace(0,udepth,201);
T_adiabat = tpotential - tabsolute + tgrad*z;

T_old = ttop + (T_adiabat(end)-ttop).*erf(z./(2*sqrt(kappa*plate_age_old)));

% Scale to adiabat
T_old = T_adiabat.*T_old./T_adiabat(end);

% Seimic oceanic lithosphere thickness at isoLAB
k = find(T_old>=isoLAB,1);
if isempty(k)
    error('Model space exceeds isoLAB')
else
    oceanic_lith_thick_old = z(k);
end

% Bottom of sediments ocean
bot_osed = sticky_layer + sed_thick;

% Upper Ocean Crust thickness (m)
upper_ocrust_thick = 3000;

% lower Ocean Crust thickness (m)
lower_ocrust_thick = 4000;

% Ocean lithosphere asthenosphere boundary (m) 
oceanLAB_old = sticky_layer + oceanic_lith_thick_old;

% CMB thermal layer depth (m)
CMB_boundarylayer = 2800e3;

% -------------------------------------------------------------------------
%% Defining lithological structure of the model
% -------------------------------------------------------------------------

% Marker counter
mm1=0;
rng('default') % seed

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
    if(MY(mm1)>bot_ocrust2 && MY(mm1)<=oceanLAB_old)
        MI(mm1) = 8;
    end
    
end

%--------------------------------------------------------------------------
% Initial temperature structure
%--------------------------------------------------------------------------

%TcmbL =  tpotential - tabsolute + tgrad*CMB_boundarylayer;  % Adiabatic temperature at CMB thermal boundary layer
TcmbL = tpotential - tabsolute;

slp_cmb = (TcmbL - tbottom)./(CMB_boundarylayer-ysize0+sticky_layer);
b_cmb = TcmbL - slp_cmb*CMB_boundarylayer;

for mm1 = 1:marknum   
    
    % Adiabatic Temperature Gradient in the Mantle
    MTK(mm1) = tpotential - tabsolute;% + tgrad*(MY(mm1) - sticky_layer);
    
    % Sticky air or sticky water
    if(MI(mm1) == 1 || MI(mm1) == 2)
        MTK(mm1) = ttop;
    end
    
    % Oceanic geotherm for subducting plate of specified age
    % after Turcotte & Schubert (2002)
    if( MY(mm1) > sticky_layer && MY(mm1)<=udepth)
       z = MY(mm1) - sticky_layer;
       m_adiabat = MTK(mm1);
       MTK(mm1) = ttop + (tpotential-tabsolute-ttop).*erf(z/(2*sqrt(kappa*plate_age_old)));
       %MTK(mm1) = ttop + (T_adiabat(end)-ttop).*erf(z/(2*sqrt(kappa*plate_age_old)));
       %MTK(mm1) = m_adiabat*MTK(mm1)/T_adiabat(end);
    end
    
    % Generate artificial thermal boundary layer near the CMB
    if MY(mm1) > CMB_boundarylayer
       MTK(mm1) = slp_cmb*MY(mm1) + b_cmb;
    end
    
    if MY(mm1) > CMB_boundarylayer - 100e3 && MX(mm1) >= 2750e3 && MX(mm1) <= 3250e3
        MTK(mm1) = MTK(mm1) + 300;
    end
    
end




