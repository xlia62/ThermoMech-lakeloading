function G = GibbsFreeEnergy(Hr,Sr,Vr,a,b,c,d,a0,k,T,P)

% Compute Gibbs Free Energy G at given Pressure (P) and Temperature (T) 
% provided with standard Pressure (1 bar) and Temperature (298K)
% Enthalpy Hr, Entropy Sr, Molar Volume Vr, 
% Cp polynomial coefficients (a,b,c,d), Thermal Expansion parameter a0 and
% Bullk Bodulus kb.  See Holland and Powell (1998)
%
% G = GibbsFreeEnergy(Hr,Sr,Vr,a,b,c,d,a0,kb,T,P)

% Standard temperature
To = 298;

% Convert Pressure to bar
P = P/1e5;

T2 = T.^2;
Ts = sqrt(T);

%  Integrate Cp and Cp/T from 298 -> T, and then mulitply by T
IntCp = a*(T-To-T.*(log(T)-log(To))) ...
    + b*(0.5*(T2 - To^2)-T.*(T-To)) ...
    - c*((1./T - 1./To)-0.5*T.*(1./T2-1./To^2)) ...
    + 2*d*((Ts-sqrt(To))+T.*(1./Ts - 1./sqrt(To)));


% Thermal expansion, 298 -> T
Vt = Vr*(1+a0*(T-To)-20*a0*(Ts-sqrt(To)));

% Compression
kt = k*(1 - 1.5e-4*(T-To));
Vp = (kt/3).*((1+4*P./kt).^(3/4) - 1);

G = Hr - T.*Sr + IntCp + Vt.*Vp;