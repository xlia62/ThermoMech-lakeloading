function [MPH410,MPH660,MI] = MantlePhaseChange(MPHASE,MTK,MPR,MI,MPH410,MPH660,MRHO,timestep,phase_function,depth)

% Function for calculating mantle phase change at 410 and 660 km depths and
% corresponding latent heat released/absorbed using univariant
% transformations (phase function is either 0 or 1)

% Phase function value for new marker tylpe
phase_change = 0.9;

% 410 km
% MPHASE(4,3) = 1.5;              % Clapeyron slope (MPa/K)
% MPHASE(5,3) = 206;              % Density change (kg/m^3)
% MPHASE(6,3) = 13.2;             % Reference phase change pressure (GPa)
% MPHASE(7,3) = 1810;             % Reference phase change temperature (K)
% MPHASE(8,3) = 1;                % Viscosity increase factor
% MPHASE(9,3) = 0.1;              % Phase transition range (0.1 GPa = 20 km);
%
% 660 km
% MPHASE(4,4) = -1.3;             % Clapeyron slope (MPa/K)
% MPHASE(5,4) = 322;              % Density change (kg/m^3)
% MPHASE(6,4) = 23.5;             % Reference phase change pressure (GPa)
% MPHASE(7,4) = 1940;             % Reference phase change temperature (K)
% MPHASE(8,4) = 10;               % Viscosity increase factor
% MPHASE(9,4) = 0.1;              % Phase transition range (0.1 GPa = 20 km);

% Caclulate the marker transition pressure using Clapeyron slope (gamma) in GPa
MP410 = MPHASE(6,3) + (MTK-MPHASE(7,3))*MPHASE(4,3)/1e3;
MP660 = MPHASE(6,4) + (MTK-MPHASE(7,4))*MPHASE(4,4)/1e3;

% Marker Presssure in GPa, 
MPR = MPR/1e9;

% locate marker types subject to phase transition
%
% Standard 14, Addons must start with 15 below
% MatProp{1} = 'Sticky Air';
% MatProp{2} = 'Sticky Water';
% MatProp{3} = 'Sediments';
% MatProp{4} = 'Upper Oceanic Crust';
% MatProp{5} = 'Lower Oceanic Crust';
% MatProp{6} = 'Upper Continental Crust';
% MatProp{7} = 'Lower Continental Crust';
% MatProp{8} = 'Oceanic Lithospheric Mantle';
% MatProp{9} = 'Continental Lithospheric Mantle';
% MatProp{10} = 'Asthenospheric Mantle';
% MatProp{11} = 'Hydrated Mantle';
% MatProp{12} = 'Upper Oceanic Plateau';
% MatProp{13} = 'Lower Oceanic Plateau';
% MatProp{14} = 'Plume Material';
% MatProp{15} = 'Eclogite';
% MatProp{16} = 'Lower Mantle'; % Perovskite
% MatProp{17} = 'PvLith';       % Equivalent to ower mantle used to track lithosphere.

k_um = MI > 7 & MI <12;
k_lm = MI > 15;
k_both = k_um | k_lm | MI == 14;

% Latent heat of phase transformation
QL410 = MPHASE(4,3)*1e6*MPHASE(5,3)*MTK./MRHO;
QL660 = MPHASE(4,4)*1e6*MPHASE(5,4)*MTK./MRHO;

% Get previous phase change fraction
oldMPH410 = MPH410(:,1);  
oldMPH660 = MPH660(:,1); 

% Calculate phase transition using eiher a univariate thransformation or a
% hyperbolic phase function

switch phase_function
    
    case 'univariate'
        
        % Univariate phase function, takes into account retrograde
        MPH410(:,1) = double(MPR >= MP410 & MPR < MP660 & k_both);  % Transizition zone
        MPH660(:,1) = double(MPR >= MP660 & k_both);                % Lower Mantle
                        
    case 'hyperbollic'
        
        % Calulate marker excess phase pressure
        MP410 = MPR - MP410;
        MP660 = MPR - MP660;
        
       
        % Caclulate phase function using the excess pressure and
        % transformation pressure range (width of transformation), it is centered
        % such that 50% is at 0 excess pressure.       
        MPH410(:,1) = 0.5*(1+tanh(MP410/MPHASE(9,3)));
        MPH660(:,1) = 0.5*(1+tanh(MP660/MPHASE(9,4)));
        
        % Ensure transformation is only for valid rock-types 
        MPH410(~k_both,1) = 0;
        MPH660(~k_both,1) = 0;                
        
end

% Caclulate transofmation latent-heat HL=QL*DG/Dt in which
% QL = gamma*drho*T/rho is the latent heat per unit volume and 
% G is the phase function for the transformation
MPH410(:,2) = (MPH410(:,1) - oldMPH410).*QL410/timestep;
MPH660(:,2) = (MPH660(:,1) - oldMPH660).*QL660/timestep;

% Transform markertype according to phase change (the depth is there to
% ensure that we don't have incorrect phase changes because of incorrect
% pressure solutions at start of simulations
k = MPH660(:,1) > phase_change & (MI== 10 | MI==11) & depth > 200e3;
MI(k) = 16;  % Perovskite

k = MPH660(:,1) > phase_change & (MI == 8 | MI== 9) & depth > 200e3;
MI(k) = 17;  % Perovskite lithosphere

% Make sure the correct prograde or retrograde transition of the 660 to upper-mantle             
k = MPH660(:,1) < phase_change & (MI == 16 | MI == 17);
MI(k) = 10;


