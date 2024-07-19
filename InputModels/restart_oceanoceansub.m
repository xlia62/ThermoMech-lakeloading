% Restart model file

%========================= Need to save this in the future runs
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
MPHASE(8,4) = 10;                % Viscosity increase factor
MPHASE(9,4) = 0.1;                % Phase transition range for hyperbolic (0.1 GPa = 20 km)
%=========================

sticky_layer = 20000;

load ../forced8/grids_27.mat
load ../forced8/markers_27.mat

startstep = 1;
timesum = 0;

tnum = length(gridt);

add_plume = false;

% Reset grid coordinates (x-dir only)
MX = MX - G.gridx(1);
G.gridx = G.gridx - G.gridx(1);

% Add ridges to left and right boundaries

% left
k = MI== 8 & MX <= G.gridx(1) + 50e3;              % Replace oceanic lithospshere portion with
MI(k) = 11;                           % Hydrated mantle weakzone

% right
k = MI== 8 & MX >= G.gridx(end) - 50e3;    % Replace oceanic lithospshere portion with
MI(k) = 11;                           % Hydrated mantle weakzone

MSXX=zeros(marknum,1);  % SIGMAxx - deviatoric normal stress, Pa
MSXY=zeros(marknum,1);  % SIGMAyy - shear stress, Pa
META=zeros(marknum,1);  % viscosity, Pa s
MEXX=zeros(marknum,1);  % EPSILONxx - normal strain rate, 1/s
MEXY=zeros(marknum,1);  % EPSILONyy - shear strain rate, 1/s
MGII=zeros(marknum,1);  % Accumulated strain
MRAT=ones(marknum,1);   % EiiMarker/EiiGrid Ratio
MXM=zeros(marknum,1);   % Cummulative Melt Fraction
MEXTC=zeros(marknum,1); % Cummulative Extracted Melt Fraction
MEXT=zeros(marknum,1);  % Extracted Melt Fraction
MLAMBDA=ones(marknum,1); % Weakening effect of melt/fluid
plastyn = zeros(marknum,1); % plastic yielding occured
