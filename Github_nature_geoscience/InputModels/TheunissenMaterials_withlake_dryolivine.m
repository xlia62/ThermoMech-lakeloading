
% -------------------------------------------------------------------------
%% Material properties
%
%  All take from Gerya (2010) except for the following rheology:
%
%  Upper and Lower Mantle (including lithosphere) taken from 
%  Garel et al, (2014), doi:10.1002/2014GC005257 and tuned to ideal 
%  viscosity profile and subduction dynamics.
%
% -------------------------------------------------------------------------

% Standard 17, Addons must start with 15 below
MatProp{1} = 'Sticky Air';
MatProp{2} = 'Sticky Water';
MatProp{3} = 'Sediments';
MatProp{4} = 'Upper Continental Crust 1'; % <---------- changed to Upper Continental Crust 1 for deformation tracking
MatProp{5} = 'Upper Continnental Crust 3'; % <---- changed to Upper Cont - measuring layer for deformation tracking
MatProp{6} = 'Upper Continental Crust 2'; % <----------- changed to Upper Contiental Crust 2 for deformation tracking
MatProp{7} = 'Lower Continental Crust';
MatProp{8} = 'Oceanic Lithospheric Mantle';
MatProp{9} = 'Continental Lithospheric Mantle';
MatProp{10} = 'Upper Mantle';
MatProp{11} = 'Weaer Upper Continetnal Crust -- Fault';   %<---- changed from hydrated mantle
MatProp{12} = 'Upper Oceanic Plateau'; % flood basalt 
MatProp{13} = 'Lower Oceanic Plateau';
MatProp{14} = 'Plume Material';
MatProp{15} = 'Eclogite';
MatProp{16} = 'Lower Mantle'; % Perovskite
MatProp{17} = 'PvLith';       % Equivalent to lower mantle used to track lithosphere.
MatProp{18} = 'Ice';       % Surface rift lake used to apply volumn load. 

% Add-ons HERE!

% MRHO  = density (kg/m3): RHO*[1-ALP*(T-273)]*[1+BET*(P-1e+5)], includes
%         melt density
% MFLOW = power-law: EPSILONii= AD*SIGMAii^n*exp[-(Ea+Va*P)/RT)
% MMU = shear modulus (Pa)
% MPL = Brittle/plastic strength (Pa): SIGMAyeild= C+sin(FI)*P (Gerya, Page 174)

%         C= C0, FI= FI0 for strain<= GAM0
%         C= C0+(C1-C0)/(GAM1-GAM0)*(strain-GAM0), FI= FI0+(FI1-FI0)/(GAM1-GAM0)*(strain-GAM0) for GAM0<strain<GAM1
%         C= C1, FI= FI1 for strain>= GAM1
% MCP = heat capacity (J/K/kg)
% MKT = thermal conductivity (W/m/K): k= k0+a/(T+77)*(1+b*P)
% MHR = radiogenic heat production (W/m^3)
% MHE = shear and adiobatic heating efficiency:  take on values 0 or 1
% MPH = phase change properties, clapyron slope, delta_rho

% Initialize MFLOW
MFLOW = zeros(18,11);

% Dummy values for diffusion creep to yield a very high viscosity in
% such that it will not contribute to composite viscosity
MFLOW(:,7) = 1;               % 1= power law (dry olivine: Kawazoe et al, 2009)
MFLOW(:,8) = 1e-42;           % AD, 1/s/MPa^n (Garelle et al 2014)
MFLOW(:,9) = 1;              % n
MFLOW(:,10) = 0;            % Ea, kJ/mol
MFLOW(:,11) = 0;             % Va, cm^3 = = = > *1e-6 for J/Pa, or *0.1 for J/bar
MFLOW(:,12) = 1;  

% Choose the experiment calibration factor for rheology, applied to all
% materials!
F = @(n) 1;  % Constant

% Axial compression
% F = @(n) 1/sqrt(3^(n+1));                      % stress-dependent viscosity, F1 in eq. 6.9 (Gerya, 2010)
% F = @(n) 1/( 2^((n-1)/n) * sqrt(3^((n+1)/n))); % strain-rate-dependent viscosity F2 in eq. 6.10 (Gerya, 2010)

% Simple shear
% F = @(n) 1/2^n;                                % stress-dependent viscosity, F1 in eq. 6.14 (Gerya, 2010)
% F = @(n) 1/( 2^((2*n-1)/n) );                  % strain-rate-dependent viscosity F2 in eq. 6.10 (Gerya, 2010)


% -------------------------------------------------------------------------
%% 1 = Weak Layer ("sticky air")
% -------------------------------------------------------------------------
MRHO(1,1) = 1.22;             % standard density, kg/m^3
MRHO(1,2) = 3e-5;             % thermal expansion, 1/K
MRHO(1,3) = 1e-11;            % compressibility, 1/Pa
MFLOW(1,1) = 0;               % 0= constant viscosity
MFLOW(1,2) = 1e+18;           % viscosity, Pa s
MFLOW(1,6) = 1;               % Lab to model factor
MMU(1,1) = 1e+20;             % shear modulus, Pa
MPL(1,1) = 0;                 % C0, Pa
MPL(1,2) = 0;                 % C1, Pa
MPL(1,3) = 0;                 % sin(FI0)
MPL(1,4) = 0;                 % sin(FI1)
MPL(1,5) = 0;                 % GAM0
MPL(1,6) = 1;                 % GAM1
MPL(1,7) = 0;                 % Maximum yield strengh
MCP(1,1) = 3000;                % Cp, J/kg
MKT(1,1) = 300;               % k0, W/m/K
MKT(1,2) = 0;                 % a, W/m
MHR(1,1) = 0;                 % radiogenic heat production, W/m^3
MHE(1,1) = 0;                 % shear and adiobatic heating efficiency

% -------------------------------------------------------------------------
%% 2 = Weak Layer ("sticky water")
% -------------------------------------------------------------------------
MRHO(2,1) = 1000;             % standard density, kg/m^3
MRHO(2,2) = 3e-5;             % thermal expansion, 1/K
MRHO(2,3) = 1e-11;            % compressibility, 1/Pa
MFLOW(2,1) = 0;               % 0= constant viscosity
MFLOW(2,2) = 1e+18;           % viscosity, Pa s
MFLOW(2,6) = 1;               % Lab to model factor
MMU(2,1) = 1e+20;               % shear modulus, Pa
MPL(2,1) = 0;                 % C0, Pa
MPL(2,2) = 0;                 % C1, Pa
MPL(2,3) = 0;                 % sin(FI0)
MPL(2,4) = 0;                 % sin(FI1)
MPL(2,5) = 0;                 % GAM0
MPL(2,6) = 1;                 % GAM1
MPL(2,7) = 0;                 % Maximum yield strengh
MCP(2,1) = 3000;                % Cp, J/kg
MKT(2,1) = 300;               % k0, W/m/K
MKT(2,2) = 0;                 % a, W/m
MHR(2,1) = 0;                 % radiogenic heat production, W/m^3
MHE(2,1) = 0;                 % shear and adiobatic heating efficiency

% -------------------------------------------------------------------------
%% 3 = Sediments
% -------------------------------------------------------------------------
MRHO(3,1) = 2600;             % standard density, GarelMaterialskg/m^3
MRHO(3,2) = 3e-5;             % thermal expansion, 1/K
MRHO(3,3) = 1e-11;            % compressibility, 1/Pa
MRHO(3,4) = 2400;             % melt density
MFLOW(3,1) = 1;               % 1= power law (wet quartzite: Ranalli, 1995)
MFLOW(3,2) = 8.574e-4;          % AD, 1/s/MPa^n
MFLOW(3,3) = 4;             % n
MFLOW(3,4) = 222.81;             % Ea, kJ/mol
MFLOW(3,5) = 0;               % Va, cm^3 = = = > *1e-6 for J/Pa, or *0.1 for J/bar
                            % uni-axial compression % Lab to model factor
n = MFLOW(3,3);                            
MFLOW(3,6) = F(n);  

MMU(3,1) = 1e+10;               % shear modulus, Pa
MPL(3,1) = 9.6e+6;              % C0, Pa
MPL(3,2) = 4e+6;              % C1, Pa
MPL(3,3) = 0.259;              % sin(FI0) strong
MPL(3,4) = 0.035;             % sin(FI1) weak
MPL(3,5) = 0;                 % GAM0
MPL(3,6) = 1;                 % GAM1
MPL(3,7) = 200e6;               % Maximum yield strength, Pa
MCP(3,1) = 1000;                % Cp, J/kg
MKT(3,1) = 0.64;              % k0, W/m/K
MKT(3,2) = 807;               % a, W/m
MHR(3,1) = 1.0e-7;            % radiogenic heat production, W/m^3
MHE(3,1) = 1;                 % shear and adiobatic heating efficiency

% -------------------------------------------------------------------------
%% 4 = Upper continental crust (granodiorite) 1
% -------------------------------------------------------------------------
MRHO(4,1) = 2750;             % standard density, kg/m^3
MRHO(4,2) = 3e-5;             % thermal expansion, 1/K
MRHO(4,3) = 1e-11;            % compressibility, 1/Pa
MRHO(4,4) = 2400;             % melt density, kg/m^3
MFLOW(4,1) = 1;               % 1= power law (wet quartzite: Ranalli, 1995)
MFLOW(4,2) = 8.574e-4;          % AD, 1/s/MPa^n
MFLOW(4,3) = 4;             % n
MFLOW(4,4) = 222.81;             % Ea, kJ/mol
MFLOW(4,5) = 0;               % Va, cm^3 = = = > *1e-6 for J/Pa, or *0.1 for J/bar
                            % uni-axial compression % Lab to model factor
n = MFLOW(4,3);                            
MFLOW(4,6) = F(n); 

MMU(4,1) = 36e+9;               % shear modulus, Pa
MPL(4,1) = 19.3e+6;              % C0, Pa,  10e6 From April
MPL(4,2) = 4e+6;              % C1, Pa,      1e6 From April
MPL(4,3) = 0.259;                 % sin(FI0)
MPL(4,4) = 0.035;                 % sin(FI1)
MPL(4,5) = 0;                 % GAM0
MPL(4,6) = 1;                 % GAM1
MPL(4,7) = 150e6;               % was 200e6 Maximum yield strength, Pa
MCP(4,1) = 1000;                % Cp, J/kg
MKT(4,1) = 0.64;              % k0, W/m/K
MKT(4,2) = 807;               % a, W/m
MHR(4,1) = 1e-7;              % was 1e-6, radiogenic heat production, W/m^3
MHE(4,1) = 1;                 % shear and adiobatic heating efficiency
% -------------------------------------------------------------------------
%% 5 = Upper continental crust 3
% -------------------------------------------------------------------------
MRHO(5,1) = 2750;             % standard density, kg/m^3
MRHO(5,2) = 3e-5;             % thermal expansion, 1/K
MRHO(5,3) = 1e-11;            % compressibility, 1/Pa
MRHO(5,4) = 2400;             % melt density, kg/m^3
MFLOW(5,1) = 1;               % 1= power law (plagioclase An75: Ranalli, 1995)
MFLOW(5,2) = 8.574e-4;          % AD, 1/s/MPa^n
MFLOW(5,3) = 4;             % n
MFLOW(5,4) = 222.81;             % Ea, kJ/mol
MFLOW(5,5) = 0;               % Va, cm^3  = = = > *1e-6 for J/Pa, or *0.1 for J/bar
                            % uni-axial compression % Lab to model factor
n = MFLOW(5,3);                            
MFLOW(5,6) = F(n);    

MMU(5,1) = 36e+9;             % shear modulus, Pa
MPL(5,1) = 19.3e+6;              % C0, Pa
MPL(5,2) = 4e+6;              % C1, Pa
MPL(5,3) = 0.259;               % sin(FI0), strong
MPL(5,4) = 0.035;               % sin(FI1), weak
MPL(5,5) = 0;                 % GAM0
MPL(5,6) = 1;                 % GAM1
MPL(5,7) = 200e6;               % Maximum yield strength, Pa
MCP(5,1) = 1000;                % Cp, J/kg
MKT(5,1) = 0.64;              % k0, W/m/K
MKT(5,2) = 807;               % a, W/m
MHR(5,1) = 1e-7;              % was 1e-6 radiogenic heat production, W/m^3
MHE(5,1) = 1;                 % shear and adiobatic heating efficiency

% -------------------------------------------------------------------------
%% 6 = Upper continental crust 2 (granodiorite)
% -------------------------------------------------------------------------
MRHO(6,1) = 2750;             % standard density, kg/m^3
MRHO(6,2) = 3e-5;             % thermal expansion, 1/K
MRHO(6,3) = 1e-11;            % compressibility, 1/Pa
MRHO(6,4) = 2400;             % melt density, kg/m^3
MFLOW(6,1) = 1;               % 1= power law (wet quartzite: Ranalli, 1995)
MFLOW(6,2) = 8.574e-4;          % AD, 1/s/MPa^n
MFLOW(6,3) = 4;             % n
MFLOW(6,4) = 222.81;             % Ea, kJ/mol
MFLOW(6,5) = 0;               % Va, cm^3 = = = > *1e-6 for J/Pa, or *0.1 for J/bar
                            % uni-axial compression % Lab to model factor
n = MFLOW(6,3);                            
MFLOW(6,6) = F(n); 

MMU(6,1) = 36e+9;               % shear modulus, Pa
MPL(6,1) = 19.3e+6;              % C0, Pa
MPL(6,2) = 4e+6;              % C1, Pa
MPL(6,3) = 0.259;       % was 0.259;                 % sin(FI0)
MPL(6,4) = 0.035;                 % sin(FI1)
MPL(6,5) = 0;                 % GAM0
MPL(6,6) = 1;                 % GAM1
MPL(6,7) = 150e6;               % was 200e6 Maximum yield strength, Pa
MCP(6,1) = 1000;                % Cp, J/kg
MKT(6,1) = 0.64;              % k0, W/m/K
MKT(6,2) = 807;               % a, W/m
MHR(6,1) = 1e-7;              %5e-7 radiogenic heat production, W/m^3
MHE(6,1) = 1;                 % shear and adiobatic heating efficiency

% -------------------------------------------------------------------------
%% 7 = Lower continental crust (diorite)
% -------------------------------------------------------------------------
MRHO(7,1) = 2900;%2900;             % standard density, kg/m^3
MRHO(7,2) = 3e-5;             % thermal expansion, 1/K
MRHO(7,3) = 1e-11;            % compressibility, 1/Pa
MRHO(7,4) = 2400;             % melt density, kg/m^3
MFLOW(7,1) = 1;               % 1= power law (plagioclase An75: Ranalli, 1995)
MFLOW(7,2) = 8.574e-4;          % AD, 1/s/MPa^n
MFLOW(7,3) = 4;             % n
MFLOW(7,4) = 222.81;             % Ea, kJ/mol
MFLOW(7,5) = 0;               % Va, cm^3 = = = > *1e-6 for J/Pa, or *0.1 for J/bar
                            % uni-axial compression % Lab to model factor
n = MFLOW(7,3);                            
MFLOW(7,6) = F(n);  
MMU(7,1) = 2.5e+10;             % shear modulus, Pa
MPL(7,1) = 19.3e+6;              % C0, Pa
MPL(7,2) = 4e+6;               % C1, Pa
MPL(7,3) = 0.259;               % sin(FI0)
MPL(7,4) = 0.035;               % sin(FI1)
MPL(7,5) = 0;                 % GAM0
MPL(7,6) = 1;                 % GAM1
MPL(7,7) = 200e6;               % Maximum yield strength, Pa
MCP(7,1) = 1000;                % Cp, J/kg
MKT(7,1) = 1.18;              % k0, W/m/K
MKT(7,2) = 474;               % a, W/m
MHR(7,1) = 1e-7;              %5.0e-7 radiogenic heat production, W/m^3
MHE(7,1) = 1;                 % shear and adiobatic heating efficiency

% -------------------------------------------------------------------------
%% 8 = Oceanic Lithospheric mantle
% -------------------------------------------------------------------------
MRHO(8,1) = 3300;             % standard density, kg/m^3
MRHO(8,2) = 3e-5;             % thermal expansion, 1/K
MRHO(8,3) = 1e-11;            % compressibility, 1/Pa

% Dislocation creep
MFLOW(8,1) = 1;               % 1= power law (dry olivine: Kawazoe et al, 2009)
MFLOW(8,2) = 1.1e+5;          % AD, 1/s/MPa^n
MFLOW(8,3) = 3.5;             % n
MFLOW(8,4) = 530;             % Ea, kJ/mol
MFLOW(8,5) = 15;              % Va, cm^3 = = = > *1e-6 for J/Pa, or *0.1 for J/bar
                              % uni-axial compression % Lab to model factor
n = MFLOW(8,3);                             
MFLOW(8,6) = F(n);

% Diffusion creep             (Same as upper mantle)
MFLOW(8,7) = 1;               % 1= power law (dry olivine: Kawazoe et al, 2009)
MFLOW(8,8) = 3e-5;           % AD, 1/s/MPa^n (Garelle et al 2014)
MFLOW(8,9) = 1;              % n
MFLOW(8,10) = 300;            % Ea, kJ/mol
MFLOW(8,11) = 4;             % Va, cm^3 = = = > *1e-6 for J/Pa, or *0.1 for J/bar
n = MFLOW(8,9);                            
MFLOW(8,12) = F(n);  

MMU(8,1) = 6.7e+10;           % shear modulus, Pa
MPL(8,1) = 2e+6;              % C0, Pa
MPL(8,2) = 2e+6;              % C1, Pa
MPL(8,3) = 0.2;               % sin(FI0)
MPL(8,4) = 0.2;              % sin(FI1)
MPL(8,5) = 0;                 % GAM0
MPL(8,6) = 1;                 % GAM1
MPL(8,7) = 200e6;             % Maximum yield strength, Pa
MCP(8,1) = 1000;              % Cp, J/kg
MKT(8,1) = 0.73;              % k0, W/m/K
MKT(8,2) = 1293;              % a, W/m
MKT(8,3) = 4e-11;             % Pressure dependance of thermal conductivity b, W/m/Pa
MHR(8,1) = 2.2e-8;              % radiogenic heat production, W/m^3
MHE(8,1) = 1;                 % shear and adiobatic heating efficiency

% -------------------------------------------------------------------------
%% 9 = Continental Lithospheric mantle -- dry mantle 
% -------------------------------------------------------------------------
MRHO(9,1) = 3300;            % standard density, kg/m^3
MRHO(9,2) = 3e-5;            % thermal expansion, 1/K
MRHO(9,3) = 1e-11;           % compressibility, 1/Pa
MRHO(9,4) = 2900;            % melt density, kg/m^3   (Gerya and Mellick, 2010)
MRHO(9,5) = -15;             % depletion (compositional density constrast) kg/m^3

MFLOW(9,1) = 1;              % 1= power law (dry olivine: Koptev et al., 2016)
% Dislocation creep
MFLOW(9,2) = 3.98e+5;    % AD, 1/s/MPa^n  3.98e16 pa^n*s (dry olivine: Koptev et al., 2016)
MFLOW(9,3) = 3.5;             % n
MFLOW(9,4) = 532;             % Ea, kJ/mol
MFLOW(9,5) = 15;              % Va, cm^3 = = = > *1e-6 for J/Pa, or *0.1 for J/bar
                              % uni-axial compression % Lab to model factor                            
n = MFLOW(9,3);                            
MFLOW(9,6) = F(n);

% Diffusion creep             (Same as upper mantle)
MFLOW(9,7) = 1;               % 1= power law (dry olivine: Kawazoe et al, 2009)
MFLOW(9,8) = 3e-45;           %was 3e-45; AD, 1/s/MPa^n (Garelle et al 2014)
MFLOW(9,9) = 1;              % n
MFLOW(9,10) = 300;            % Ea, kJ/mol
MFLOW(9,11) = 4;             % Va, cm^3 = = = > *1e-6 for J/Pa, or *0.1 for J/bar
n = MFLOW(9,9);                            
MFLOW(9,12) = F(n);  

MMU(9) = 74e+9;              % shear modulus, Pa
MPL(9,1) = 19.3e+6;              % C0, Pa
MPL(9,2) = 4e+6;              % C1, Pa
MPL(9,3) = 0.259;                % sin(FI0)
MPL(9,4) = 0.035;                % sin(FI1)
MPL(9,5) = 0;                  % GAM0
MPL(9,6) = 1;                  % GAM1
MPL(9,7) = 200e6;               % Maximum yield strength, Pa
MCP(9,1) = 1000;               % Cp, J/kg
MKT(9,1) = 0.73;               % k0, W/m/K  (Gerya and Mellick, 2010)
MKT(9,2) = 1293;               % a, W/m
MKT(9,3) = 4e-11;             % Pressure dependance of thermal conductivity b, W/m/Pa
MHR(9,1) = 4e-9;             % radiogenic heat production, W/m^3
MHE(9,1) = 1;                  % shear and adiobatic heating efficiency

% -------------------------------------------------------------------------
%% 10 = Upper mantle -- dry mantle 2
% -------------------------------------------------------------------------

MRHO(10,1) = 3300;            % standard density, kg/m^3
MRHO(10,2) = 3e-5;            % thermal expansion, 1/K
MRHO(10,3) = 1e-11;           % compressibility, 1/Pa
MRHO(10,4) = 2900;            % melt density, kg/m^3   (Gerya and Mellick, 2010)
MRHO(10,5) = -15;             % depletion (compositional density constrast) kg/m^3

MFLOW(10,1) = 1;              % 1= power law (dry olivine: Koptev et al., 2016)
% Dislocation creep
MFLOW(10,2) = 3.98e+5;    % AD, 1/s/MPa^n  3.98e16 pa^n*s (dry olivine: Koptev et al., 2016)
MFLOW(10,3) = 3.5;             % n
MFLOW(10,4) = 532;             % Ea, kJ/mol
MFLOW(10,5) = 15;              % Va, cm^3 = = = > *1e-6 for J/Pa, or *0.1 for J/bar
                              % uni-axial compression % Lab to model factor                            
n = MFLOW(10,3);                            
MFLOW(10,6) = 0.3;            % change visocity to 1e19 - 1e20 level

% Diffusion creep             (Same as upper mantle)
MFLOW(10,7) = 1;               % 1= power law (dry olivine: Kawazoe et al, 2009)
MFLOW(10,8) = 3e-45;           %was 3e-45; AD, 1/s/MPa^n (Garelle et al 2014)
MFLOW(10,9) = 1;              % n
MFLOW(10,10) = 300;            % Ea, kJ/mol
MFLOW(10,11) = 4;             % Va, cm^3 = = = > *1e-6 for J/Pa, or *0.1 for J/bar
n = MFLOW(10,9);                            
MFLOW(10,12) = 0.3;  

MMU(10) = 74e+9;              % shear modulus, Pa
MPL(10,1) = 19.3e+6;              % C0, Pa
MPL(10,2) = 4e+6;              % C1, Pa
MPL(10,3) = 0.259;                % sin(FI0)
MPL(10,4) = 0.035;                % sin(FI1)
MPL(10,5) = 0;                  % GAM0
MPL(10,6) = 1;                  % GAM1
MPL(10,7) = 200e6;               % Maximum yield strength, Pa
MCP(10,1) = 1000;               % Cp, J/kg
MKT(10,1) = 0.73;               % k0, W/m/K  (Gerya and Mellick, 2010)
MKT(10,2) = 1293;               % a, W/m
MKT(10,3) = 4e-11;             % Pressure dependance of thermal conductivity b, W/m/Pa
MHR(10,1) = 4e-9;             % radiogenic heat production, W/m^3
MHE(10,1) = 1;                  % shear and adiobatic heating efficiency

% -------------------------------------------------------------------------
%% 11 = Upper continental crust (granodiorite) -- weak faulting
% -------------------------------------------------------------------------
MRHO(11,1) = 2750;             % standard density, kg/m^3
MRHO(11,2) = 3e-5;             % thermal expansion, 1/K
MRHO(11,3) = 1e-11;            % compressibility, 1/Pa
MRHO(11,4) = 2400;             % melt density, kg/m^3
MFLOW(11,1) = 1;               % 1= power law (wet quartzite: Ranalli, 1995)
MFLOW(11,2) = 8.574e-4;          % AD, 1/s/MPa^n
MFLOW(11,3) = 4;             % n
MFLOW(11,4) = 222.81;             % Ea, kJ/mol
MFLOW(11,5) = 0;               % Va, cm^3 = = = > *1e-6 for J/Pa, or *0.1 for J/bar
                            % uni-axial compression % Lab to model factor
n = MFLOW(11,3);                            
MFLOW(11,6) = F(n); 

MMU(11,1) = 36e+9;               % shear modulus, Pa

% MPL = Brittle/plastic strength (Pa), major difference from 6: SIGMAyeild= C+sin(FI)*P
%         C= C0, FI= FI0 for strain<= GAM0
%         C= C0+(C1-C0)/(GAM1-GAM0)*(strain-GAM0), FI= FI0+(FI1-FI0)/(GAM1-GAM0)*(strain-GAM0) for GAM0<strain<GAM1
%         C= C1, FI= FI1 for strain>= GAM1

MPL(11,1) = 4e+6;              % was 4e6, C0, Pa, difference
MPL(11,2) = 4e+6;              % was 4e6, C1, Pa
MPL(11,3) = 0.035;                 % sin(FI0) difference
MPL(11,4) = 0.035;                 % sin(FI1)
MPL(11,5) = 0;                   % GAM0
MPL(11,6) = 1;                   % GAM1
MPL(11,7) = 200e6;               % Maximum yield strength, Pa
MCP(11,1) = 1000;                % Cp, J/kg
MKT(11,1) = 0.64;              % k0, W/m/K
MKT(11,2) = 807;               % a, W/m
MHR(11,1) = 1.0e-7;              % was 1e-6, radiogenic heat production, W/m^3
MHE(11,1) = 1;                 % shear and adiobatic heating efficiency


% -------------------------------------------------------------------------
%% 12 = Upper oceanic crust plateau (basalts)
% -------------------------------------------------------------------------
MRHO(12,1) = 3000;         % standard density, kg/m^3
MRHO(12,2) = 3e-5;             % thermal expansion, 1/K
MRHO(12,3) = 1e-11;            % compressibility, 1/Pa
MRHO(12,4) = 2400;             % melt density kg/m^3
MFLOW(12,1) = 1;               % 1= power law (wet quartzite: Ranalli, 1995)
MFLOW(12,2) = 3.2e-4;          % AD, 1/s/MPa^n
MFLOW(12,3) = 2.3;             % n
MFLOW(12,4) = 154;             % Ea, kJ/mol
MFLOW(12,5) = 0;               % Va, cm^3 = = = > *1e-6 for J/Pa, or *0.1 for J/bar
                            % uni-axial compression % Lab to model factor
n = MFLOW(12,3);                            
MFLOW(12,6) = F(n);  

MMU(12,1) = 2.5e+10;             % shear modulus, Pa
MPL(12,1) = 1e+6;              % C0, Pa
MPL(12,2) = 1e+6;              % C1, Pa
MPL(12,3) = 0.03;                 % sin(FI0)
MPL(12,4) = 0.03;                 % sin(FI1)
MPL(12,5) = 0;                   % GAM0
MPL(12,6) = 0.1;                 % GAM1
MPL(12,7) = 200e6;               % Maximum yield strength, Pa
MCP(12,1) = 1000;                % Cp, J/kg
MKT(12,1) = 1.18;              % k0, W/m/K
MKT(12,2) = 474;               % a, W/m
MHR(12,1) = 2.5e-7;              % radiogenic heat production, W/m^3
MHE(12,1) = 1;                 % shear and adiobatic heating efficiency

% -------------------------------------------------------------------------
%% 13 = Lower oceanic crust plateau
% -------------------------------------------------------------------------
MRHO(13,1) = 3000-150;             % standard density, kg/m^3
MRHO(13,2) = 3e-5;             % thermal expansion, 1/K
MRHO(13,3) = 1e-11;            % compressibility, 1/Pa
MFLOW(13,1) = 1;               % 1= power law (plagioclase An75: Ranalli, 1995)
MFLOW(13,2) = 3.3e-4;          % AD, 1/s/MPa^n
MFLOW(13,3) = 3.2;             % n
MFLOW(13,4) = 238;             % Ea, kJ/mol
MFLOW(13,5) = 0;               % Va, cm^3 = = = > *1e-6 for J/Pa, or *0.1 for J/bar
                            % uni-axial compression % Lab to model factor
n = MFLOW(13,3);                            
MFLOW(13,6) = F(n);  

MMU(13,1) = 2.5e+10;             % shear modulus, Pa
MPL(13,1) = 1e+6;              % C0, Pa
MPL(13,2) = 1e+6;              % C1, Pa
MPL(13,3) = 0.2;               % sin(FI0)
MPL(13,4) = 0.2;               % sin(FI1)
MPL(13,5) = 0;                 % GAM0
MPL(13,6) = 1;                 % GAM1
MPL(13,7) = 200e6;               % Maximum yield strength, Pa
MCP(13,1) = 1000;                % Cp, J/kg
MKT(13,1) = 1.18;              % k0, W/m/K
MKT(13,2) = 474;               % a, W/m
MHR(13,1) = 2.5e-7;              % radiogenic heat production, W/m^3
MHE(13,1) = 1;                 % shear and adiobatic heating efficiency

% -------------------------------------------------------------------------
%% 14 Plume material tracking equivlent to Asthenosphere = 6
% -------------------------------------------------------------------------

MRHO(14,1) = 3300;             % standard density, kg/m^3
MRHO(14,2) = 3e-5;             % thermal expansion, 1/K
MRHO(14,3) = 1e-11;            % compressibility, 1/Pa
MRHO(14,4) = 2900;             % melt density, kg/m^3
MRHO(14,5) = 10;              % depleted

% Disclocation creep
MFLOW(14,1) = 1;               % 1= power law (dry olivine: Kawazoe et al, 2009)
MFLOW(14,2) = 1.1e+5;          % AD, 1/s/MPa^n
MFLOW(14,3) = 3.5;             % n
MFLOW(14,4) = 530;             % Ea, kJ/mol
MFLOW(14,5) = 0;              % Va, cm^3 = = = > *1e-6 for J/Pa, or *0.1 for J/bar
                              % uni-axial compression % Lab to model factor                            
n = MFLOW(10,3);                            
MFLOW(14,6) = F(n);  

% Diffusion creep        
MFLOW(14,7) = 1;               % 1= power law (dry olivine: Kawazoe et al, 2009)
MFLOW(14,8) = 3e-5;           % AD, 1/s/MPa^n (Garelle et al 2014)
MFLOW(14,9) = 1;              % n
MFLOW(14,10) = 300;            % Ea, kJ/mol
MFLOW(14,11) = 4;             % Va, cm^3 = = = > *1e-6 for J/Pa, or *0.1 for J/bar
n = MFLOW(14,9);                            
MFLOW(14,12) = F(n);  

MMU(14,1) = 6.7e+10;           % shear modulus, Pa
MPL(14,1) = 2e+6;              % C0, Pa
MPL(14,2) = 2e+6;              % C1, Pa
MPL(14,3) = 0.2;               % sin(FI0) friction coefficient
MPL(14,4) = 0.2;               % sin(FI1)
MPL(14,5) = 0;                 % GAM0
MPL(14,6) = 1;                 % GAM1
MPL(14,7) = 200e6;               % Maximum yield strength, Pa
MCP(14,1) = 1000;                % Cp, J/kg
MKT(14,1) = 0.73;              % k0, W/m/K
MKT(14,2) = 1293;              % a, W/m
MKT(14,3) = 4e-11;             % Pressure dependance of thermal conductivity b, W/m/Pa
MHR(14,1) = 2.2e-8;              % radiogenic heat production, W/m^3
MHE(14,1) = 1;                 % shear and adiobatic heating efficiency

% -------------------------------------------------------------------------
%% 15 = Eclogite (Use Gabbro Rheology instead of basalt)
% -------------------------------------------------------------------------
MRHO(15,1) = 3400;             % standard density, kg/m^3
MRHO(15,2) = 3e-5;             % thermal expansion, 1/K
MRHO(15,3) = 1e-11;            % compressibility, 1/Pa
MRHO(15,4) = 2900;             % melt density, kg/m^3
MFLOW(15,1) = 1;               % 1= power law (plagioclase An75: Ranalli, 1995)
MFLOW(15,2) = 3.3e-4;          % AD, 1/s/MPa^n
MFLOW(15,3) = 3.2;             % n
MFLOW(15,4) = 480;             % Ea, kJ/mol
MFLOW(15,5) = 0;               % Va, cm^3 = = = > *1e-6 for J/Pa, or *0.1 for J/bar
                            % uni-axial compression % Lab to model factor
n = MFLOW(15,3);                            
MFLOW(15,6) = F(n);  

MMU(15,1) = 2.5e+10;           % shear modulus, Pa
MPL(15,1) = 1e+6;              % C0, Pa
MPL(15,2) = 1e+6;              % C1, Pa
MPL(15,3) = 0.6;               % sin(FI0), strong
MPL(15,4) = 0.3;               % sin(FI1), weak
MPL(15,5) = 0;                 % GAM0
MPL(15,6) = 1;                 % GAM1
MPL(15,7) = 200e6;               % Maximum yield strength, Pa
MCP(15,1) = 1000;              % Cp, J/kg
MKT(15,1) = 1.18;              % k0, W/m/K
MKT(15,2) = 474;               % a, W/m
MHR(15,1) = 2.5e-7;            % radiogenic heat production, W/m^3
MHE(15,1) = 1;                 % shear and adiobatic heating efficiency

% -------------------------------------------------------------------------
%% 16 = Perovskite 
%  Keep density same as upper mantle, only change rheology because the
%  density increase is handled through the phase function in the equation
%  of state
% -------------------------------------------------------------------------
MRHO(16,:) = MRHO(10,:);

% Dislocation creep D
% use dummy AD so it does not occure in lower mantle
MFLOW(16,1) = 1;            % 1= power law ( after Garell et al, 2014)
MFLOW(16,2) = 1e-21;        % AD, 1/s/MPa^n to convert from Pa^n, (1e-6*AD^(-1/n))^-n
MFLOW(16,3) = 3.5;             % n
MFLOW(16,4) = 530;             % Ea, kJ/mol
MFLOW(16,5) = 15;              % Va, cm^3 = = = > *1e-6 for J/Pa, or *0.1 for J/bar
                              % uni-axial compression % Lab to model factor  
n = MFLOW(16,3);
MFLOW(16,6) = F(n);            % uni-axial compression % Lab to model factor

% Diffusion creep
MFLOW(16,7) = 1;            % 1= power law ( after Garell et al, 2014)
MFLOW(16,8) = 2e-10;      % AD, 1/s/MPa^n to convert from Pa^n, (1e-6*AD^(-1/n))^-n
MFLOW(16,9) = 1;            % n
MFLOW(16,10) = 200;          % Ea, kJ/mol
MFLOW(16,11) = 1.5;          % Va, cm^3 = = = > *1e-6 for J/Pa, or *0.1 for J/bar
n = MFLOW(16,9);
MFLOW(16,12) = F(n);             % uni-axial compression % Lab to model factor

MMU(16,1) = 6.7e+10;          % shear modulus, Pa
MPL(16,1) = 2e+6;            % C0, Pa
MPL(16,2) = 2e+6;             % C1, Pa
MPL(16,3) = 0.2;               % sin(FI0)
MPL(16,4) = 0.2;               % sin(FI1)
MPL(16,5) = 0;                 % GAM0
MPL(16,6) = 1;                 % GAM1
MPL(16,7) = 100e6;             % Maximum yield strength, Pa

MCP(16,:) = MCP(10,:);
MKT(16,:) = MKT(10,:);
MHR(16,:) = MHR(10,:);
MHE(16,:) = MHE(10,:);

% -------------------------------------------------------------------------
%% 17 = Perovskite Lithosphere 
% Only used for tracking, keep same oceanic lithosphere and change rheology
% to perovskite (type 16)
% -------------------------------------------------------------------------
MRHO(17,:) = MRHO(8,:);
MFLOW(17,:) = MFLOW(16,:);
MMU(17,:) = MMU(8,:);
MPL(17,:) = MPL(8,:);
MCP(17,:) = MCP(8,:);
MKT(17,:) = MKT(8,:);
MHR(17,:) = MHR(8,:);
MHE(17,:) =  MHE(8,:);

% OVERWRITE YIELD STRENGHT (Backward Compatibility)
MPL(:,7) = maxyield;

% -------------------------------------------------------------------------
%% 18 = Rift Lake 
% -------------------------------------------------------------------------

% shear and adiobatic heating efficiency
MRHO(18,1) = 1000;             % standard density, kg/m^3
MRHO(18,2) = 3e-5;             % thermal expansion, 1/K
MRHO(18,3) = 1e-11;            % compressibility, 1/Pa
% MFLOW(18,1) = 1;               % 1= power law (plagioclase An75: Ranalli, 1995) EPSILONii= AD*SIGMAii^n*exp[-(Ea+Va*P)/RT)
% MFLOW(18,2) = 1e-4;          % AD, 1/s/MPa^n
% MFLOW(18,3) = 4;             % n
% MFLOW(18,4) = 222.81;         % Ea, kJ/mol
% MFLOW(18,5) = 0;               % Va, cm^3 = = = > *1e-6 for J/Pa, or *0.1 for J/ba
% n = MFLOW(18,3);                            
% MFLOW(18,6) = F(n); 
MFLOW(18,1) = 0;               % 0= constant viscosity
MFLOW(18,2) = 1e+18;           % viscosity, Pa s
MFLOW(18,6) = 1;               % Lab to model factor
MMU(18,1) = 1e+20;             % shear modulus, Pa, Titanium-like
MPL(18,1) = 0;              % C0, Pa
MPL(18,2) = 0;               % C1, Pa
MPL(18,3) = 0;               % sin(FI0)
MPL(18,4) = 0;               % sin(FI1)
MPL(18,5) = 0;                 % GAM0
MPL(18,6) = 1;                 % GAM1
MPL(18,7) = 0;               % Maximum yield strength, Pa
MCP(18,1) = 3200;                % Cp, J/kg
MKT(18,1) = 300;              % k0, W/m/K
MKT(18,2) = 0;               % a, W/m
MHR(18,1) = 0;              % radiogenic heat production, W/m^3
MHE(18,1) = 0;                 % shear and adiobatic heating efficiency


