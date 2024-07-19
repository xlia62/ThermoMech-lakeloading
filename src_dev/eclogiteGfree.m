function dG = eclogiteGfree(T,P)

% Calculate the change in Gibbs free energy for the following solid state
% reaction:
%
%           Albite ==> Jadite + Quartz
%

% Albite NaAlSi3O8  (1 bar, T = 25oC)
dHr = -3934600;         % enthalpy J
Sr = 210.1;             % entropy J
Vmol = 10.006;          % Molar volume J/bar
a = 452;        % Heat Capicity Term 1, J/K
b = -0.013364;  % Heat Capicity Term 2, J
c = -1275900;   % Heat Capicity Term 3, J.K
d = -3953.6;    % Heat Capicity Term 4, J/(K^1/2)
alpha0 = 0.0000445;     % Thermal Expansion
kbulk = 593000;         % Bulk Modulus

G1 = GibbsFreeEnergy(dHr,Sr,Vmol,a,b,c,d,alpha0,kbulk,T,P);

% Jadite NaAlSi2O6  (1 bar, T = 25oC)
dHr = -3027830;         % enthalpy J
Sr = 133.5;             % entropy J
Vmol = 6.04;            % Molar volume J/bar
a = 301.1;      % Heat Capicity Term 1, J/K
b = 0.010143;  % Heat Capicity Term 2, J
c = -2239300;   % Heat Capicity Term 3, J.K
d = -2055.1;    % Heat Capicity Term 4, J/(K^1/2)
alpha0 = 0.0000466;     % Thermal Expansion
kbulk = 1284000;        % Bulk Modulus

G2 = GibbsFreeEnergy(dHr,Sr,Vmol,a,b,c,d,alpha0,kbulk,T,P);

% Quartz SiO2 (1 bar, T = 25oC)
dHr = -910880;          % enthalpy J
Sr = 41.5;              % entropy J
Vmol = 2.269;           % Molar volume J/bar
a = 110.7;      % Heat Capicity Term 1, J/K
b = -0.005189;  % Heat Capicity Term 2, J
c = 9;          % Heat Capicity Term 3, J.K
d = -1128.3;    % Heat Capicity Term 4, J/(K^1/2)
alpha0 = 0.0000065;     % Thermal Expansion
kbulk = 750000;         % Bulk Modulus

G3 = GibbsFreeEnergy(dHr,Sr,Vmol,a,b,c,d,alpha0,kbulk,T,P);

% Albeit ==> Jadite + Quartz
% Products - Reactants
dG = (G2+G3) - G1;

