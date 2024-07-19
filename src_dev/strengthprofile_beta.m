function S = strengthprofile(t1,rho1,flow1,pla1,...
    t2,rho2,flow2,pla2,...
    t3,rho3,flow3,pla3,...
    t4,rho4,flow4,pla4,...
    z5,rho5,flow5,pla5,S)

% Generate a strength profile for an upper 5-layer Earth, layer 1 is
% sediments, layer 2 is upper-crust, layer 3 is lower-crust, layer 4 is
% mantle lithosphere, and layer 5 is the asthenosphere
%
% S = strengthprofile(t1,rho1,flow1,pla1,...
%                     t2,rho2,flow2,pla2,...
%                     t3,rho3,flow3,pla3,...
%                     t4,rho4,flow4,pla4,...
%                     z5,rho5,flow5,pla5,...
%                     S)
% 
% where:
%      t1,...,t4 correspond to thickness of each layer 1,...,4 and z5 is
%      the total depth of the model, t1+t2+t3+t4 <= z5
%      rho1,...,rho5 density vector for each layer in which
%           rho(1) = reference density
%           rho(2) = thermal expansivity
%           rho(3:4) = not used, set to zero
%           rho(5) = chemical density component (negative for depletion)
%      flow1,...,flow5 rheology vector for each layer in which
%           flow(1) = flow law, 0 for constant, 1 for power loaw
%           flow(2) = Ad coefficient for power law (1/s/MPa^n), viscosity
%                     for constant (Pa.s)
%           flow(3) = power law exponent
%           flow(4) = Activation Energy kJ/mol
%           flow(5) = Volume Activation Energy cm^3
%      pla1,....,pla5 brittle Draker-Prager constants for each layer
%           pla(1) = C0 cohesion (Pa)
%           pla(3) = friction coefficient, mu or tan(FI0)
%           pla(7) = von Mises yield criterion (Pa), if pla(7) does not
%                    exist, will use inf
%           other entries in pla vector are ignored, set to zero

display ('test');
% gravity
g = 9.81;
R = 8.314;    % Gas constant
eii = S.eii;

% 1st layer; sediments
if t1 > 0
    k = S.z <= t1;
    z = S.z(k);
    T = S.T(k);
    bnd1 = z(end);
    
    % Equation of state (assume temperature is dominating effect)
    den1 = rho1(1)*(1-rho1(2)*(T-273)) + rho1(5);  
    
    o1 = ones(size(den1));
 
else
    den1 = [];
    bnd1 = [];
    o1 = [];
    
end

% Rheological factors
n1 = flow1(3)*o1;
Ad1 = flow1(2)*o1;
Ea1 = flow1(4)*o1;
Va1 = flow1(5)*o1;
C01 = pla1(1)*o1;
F01 = pla1(3)*o1;
if length(pla1)>6
    M01 = pla1(7)*o1;
else
    M01 = inf*o1;
end

% 2nd layer, upper crust
k = t1 < S.z & S.z <= t1+t2;
z = S.z(k);
T = S.T(k);
bnd2 = z(end);
% Equation of state (assume temperature is dominating effect)
den2 = rho2(1)*(1-rho2(2)*(T-273))+ rho2(5);

% Rheological factors
o2 = ones(size(den2));
n2 = flow2(3)*o2;
Ad2 = flow2(2)*o2;
Ea2 = flow2(4)*o2;
Va2 = flow2(5)*o2;
C02 = pla2(1)*o2;
F02 = pla2(3)*o2;
if length(pla2)>6
    M02 = pla2(7)*o2;
else
    M02 = inf*o2;
end

% 3rd layer, lower crust
k = t1+t2 < S.z & S.z <= t3+t2+t1;
z = S.z(k);
T = S.T(k);
bnd3 = z(end);
% Equation of state (assume temperature is dominating effect)
den3 = rho3(1)*(1-rho3(2)*(T-273))+ rho3(5);

% Rheological factors
o3 = ones(size(den3));
n3 = flow3(3)*o3;
Ad3 = flow3(2)*o3;
Ea3 = flow3(4)*o3;
Va3 = flow3(5)*o3;
C03 = pla3(1)*o3;
F03 = pla3(3)*o3;
if length(pla3)>6
    M03 = pla3(7)*o3;
else
    M03 = inf*o3;
end

% 4th layer, mantle lithosphere
k = t1+t2+t3 < S.z & S.z <= t4+t3+t2+t1;
z = S.z(k);
T = S.T(k);
bnd4 = z(end);
% Equation of state (assume temperature is dominating effect)
den4 = rho4(1)*(1-rho4(2)*(T-273))+ rho4(5);

% Rheological factors
o4 = ones(size(den4));
Ad4 = flow4(2)*o4;
n4 = flow4(3)*o4;
Ea4 = flow4(4)*o4;
Va4 = flow4(5)*o4;
C04 = pla4(1)*o4;
F04 = pla4(3)*o4;
if length(pla4)>6
    M04 = pla4(7)*o4;
else
    M04 = inf*o4;
end

% 5th layer, asthenospheric mantle 
k = t1+t2+t3+t4 < S.z & S.z <= z5;
z = S.z(k);
if isempty(z)
    error('Increase depth of model space, z5')
end
T = S.T(k);
bnd5 = z(end);

% Equation of state (assume temperature is dominating effect)
den5 = rho5(1)*(1-rho5(2)*(T-273))+ rho5(5);

% Rheological factors
o5 = ones(size(den5));
Ad5 = flow5(2)*o5;
n5 = flow5(3)*o5;
Ea5 = flow5(4)*o5;
Va5 = flow5(5)*o5;
C05 = pla5(1)*o5;
F05 = pla5(3)*o5;
if length(pla5)>6
    M05 = pla5(7)*o5;
else
    M05 = inf*o5;
end

bnd = [bnd1;bnd2;bnd3;bnd4;bnd5];
rho = [den1;den2;den3;den4;den5];
n = [n1;n2;n3;n4;n5];
Ad = [Ad1;Ad2;Ad3;Ad4;Ad5];
Ea = [Ea1;Ea2;Ea3;Ea4;Ea5];
Va = [Va1;Va2;Va3;Va4;Va5];
C0 = [C01;C02;C03;C04;C05];
F0 = [F01;F02;F03;F04;F05];
M0 = [M01;M02;M03;M04;M05];

% Rheological boundaries
S.bounds = bnd;
% Caclulate pressure
S.P = g*cumtrapz(S.z,rho);

% Lithospheric strength, 2nd invarient stress
% Compute the exponential term
syield = C0 + F0.*S.P;  % Pa
sii = ((eii./Ad).^(1./n)).*exp((Ea*1000.+Va.*S.P*1e-6)./(n*R.*S.T)); % MPa

% The case of constant viscosity
k = n == 0;
sii(k) = 2*eii*Ad(k);

S.sigd = min(sii*1e6,syield); % Pa
S.sigd = min(S.sigd,M0);  % von Mises criterion
