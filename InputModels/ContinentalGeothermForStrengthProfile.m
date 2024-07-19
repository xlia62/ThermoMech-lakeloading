%This function can be used to measure the temperature at any depth Z.

function  T = ContinentalGeothermForStrengthProfile(n,Z,k,hr,S0, qm, Zrad, Ti, Tl, Zl,tgrad,Tp)
% n = 0 => piece-wise linear
% n = 1 => exponential decrease of heat production
% n = 2 => Constant heat production
% Need to add a equation to change the surface tempearture. DEfaults value
% is 0 Celcius

% k = 3.5; %k - Thermal Conductivity (not used for n = 0)
% hr = 10000; % hr - Characteristic drop off for n > 1, Moho depth for n=0
% S0 = 3e-6; % S0 - Surface heat production rate for n > 1, Moho
%                   temperature for n = 0
% qm = 0.03;
% Zrad = 7000;% Constant Radioactive zone  (not used for n=0 or n=1)
% Ti = 0; %Ti - Surface Temperature in Celcius
% Tl = 1300; % Tl = Temperature of the base of Mantle Lithosphere
% Zl = 100000; % Zl = Depth of the base of the mantle lithosphere
%
% tgrad is the adiabatic gradient in oC/m
% Tp = mantle potential temperature
% Temperature

T = zeros(size(Z));
N = length(T);

switch n
    
    case 0 % Calculate continental geotherm using a piece-wise linear approximation
        
        % Define moho temperature
        Tmoho = S0;
        
        % Define moho depth
        moho = hr;
        
        % Linear Tempearture in crust
        i = Z <= moho;
        
        % Slope of geotherm
        m = (Tmoho-Ti)/moho;
        b = Tmoho - m*moho;
        
        T(i) = m*Z(i) + b;
        
        % Moho to base of lithosphere, intersection with adiabat
        % Temperature of adiabat at desired depth
        %Tl = Tp + tgrad*Zl;

        % Find points in the mantle lithosphere
        i = Z > moho; % & Z <= Zl;
        
        % Slope of geotherm
        m = (Tl-Tmoho)/(Zl-moho);
        b = Tl - m*Zl;
        
        T(i) = m*Z(i) + b;
        
    case 1
        
        % Lithosphere exponential decay
        i = Z <= Zl;
        T(i) = Ti + Z(i)*(Tl-Ti)/Zl + (hr^2*S0/k)*(1 - exp(-Z(i)/hr) - (Z(i)/Zl)*(1 - exp(-Zl/hr)));       

    case 2
        
        % Radioactive zone
        i = Z <= Zrad;
        T(i) = Ti + (S0*Z(i).*(Zrad - Z(i)/2) + qm*Z(i))/k;                 

        i = Z > Zrad & Z <= Zl;
        T(i) = Ti + (S0*(Zrad^2)/2 + qm*Z(i))/k;
        
end

if N > 1
    
        % Extrapolate temperature across entire domain to intersect adiabat
        i = Z <= Zl;
        j = Z > Zl;
        
        t = T(i);
        z = Z(i);
        % Caclulate the slope 
        if sum(j) > 0
            m = (t(end)-t(end-1))/(z(end)-z(end-1));
            b = t(end) - m*z(end);
            T(j) = m*Z(j) + b;
        end

end

% Now we limit the temperature it to the adiabatic gradient
Tadiabat = Tp + tgrad*Z;
j = T > Tadiabat;
T(j) = Tadiabat(j);

% j = Z >= Zl
% T(j) = Tadiabat(j);