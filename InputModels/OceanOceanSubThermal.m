function Ti = OceanOceanSubThermal(Xi,Zi,age_sp,age_op,x1,x2,thetas,Tpotential,rcurve,tgrad)

% Plate speed of age_sp
u1 = x1/(age_sp*1e6); % km/yr

age1 = 0:.1:age_sp;  % Myr
age1 = age1*1e6;       % yr

x1 = u1*age1;        % km
nx = length(x1);

z = linspace(0,rcurve,201);
nz = length(z);

T1 = zeros(nz,nx);
for i = 1:nx
    T1(:,i) = HeatEq1D(0,Tpotential,z*1000,age1(i),tgrad);
end

u2 = x2/(age_op*1e6);

age2 = 0.1:0.1:age_op;
age2 = age2*1e6;

x2 = u2*age2;
nx2 = length(x2);

T2 = zeros(nz,nx2);
for i = 1:nx2
    T2(:,i) = HeatEq1D(0,Tpotential,z*1000,age2(i),tgrad);
end

x = [x1 x2 + x1(end)];
T = [T1 fliplr(T2)];
[X,Z]=meshgrid(x,z);

ntheta = length(thetas);
b = [zeros(size(z));z-z(end)];

Xp = zeros(length(z),ntheta);
Zp = Xp;
Tp = Xp;

for i = 1:ntheta
    
    theta = thetas(i);
    R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
    bp = R*b;
    xp = bp(1,:)+x1(end);
    zp = bp(2,:)+z(end);
    
    Xp(:,i) = xp(:);
    Zp(:,i) = zp(:);
    Tp(:,i) = T1(:,end);
end

F = scatteredInterpolant(Xp(:),Zp(:),Tp(:),'linear','none');
Ts = F(X,Z);
k = isnan(Ts);
Ts(k) = T(k);

xx = linspace(0,max(x),1000);
[XX,ZZ]=meshgrid(xx,z);
TT = interp2(X,Z,Ts,XX,ZZ);
TTs = smooth2a(TT,3,3);
% for i = 1:2
% TTs = smooth2a(TTs,3,3);
% end
Ti = interp2(XX,ZZ,TTs,Xi,Zi);