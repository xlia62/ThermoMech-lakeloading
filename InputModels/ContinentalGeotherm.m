%This function can be used to measure the temperature at any Z Z.

function  T = ContinentalGeotherm(n,Z,k,hr,S0, qm, Zrad, Ti, Tl, Zl)
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
% Tl = 1200; % Tl = Temperature of the base of Mantle Lithosphere
% Zl = 100000; % Zl = Depth of the base of the mantle lithosphere


switch n
    
    case 0 % Calculate continental geotherm using a piece-wise linear approximation
        
        % Define moho temperature
        Tmoho = S0;
        
        % Define moho depth
        moho = hr;
        
        % Temperature
        T = zeros(size(Z));
        
        % Linear Tempearture in crust
        k = Z <= moho;
        nk = sum(k);
        if nk == 0
            error('Increase the depth of Moho')
        end
        T(k) = linspace(Ti,Tmoho,nk);
        
        % Moho to lithosphere, to adiabat
        k = Z >= moho & Z <= z1;
        nk = sum(k);
        if nk == 0
            error('Increase the lithospheric depth')
        end
        tmp = linspace(Tmoho,Tl,nk);
        z = Z(k);
        T(k) = tmp;
        
        k = Z > Zl;
        if sum(k) > 0 % Extrapolate temperature across entire domain
            m = (tmp(2) - tmp(1))/(z(2)-z(1));
            b = Tl - m*z(end);
            T(k) = m*Z(k) + b;
        end
        
        % Now let's bind this geotherm by the adiabat
        Tm = Ti + Tpotential + Tgrad*Z/1000;
        
        k = T > Tm;
        T(k) = Tm(k);
        
    case 1  % Constant heat production
        
        T = (Z*Tl)/Zl+(((hr^2)*S0)/k)*((1-exp(-(Z/hr)))-(Z/Zl)*(1-exp(-(Zl/hr))))+(Ti*(1-(Z/Zl)));
        
    case 2  % Exponential heat production
        
        if (Z<=Zrad)
            T = Ti+(S0*Z*(Zrad-(Z/2))+(qm*Z))/k;
        else
            T = Ti+(((S0*(Zrad^2))/2)+(qm*Zrad)+qm*(Z-Zrad))/k;
        end
end

end
