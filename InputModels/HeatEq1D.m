
function T = HeatEq1D(To,Tp,z,age,gradT)

% Thermal diffusivity
kappa = 1;

% Convert kappa into m^2/yr
kappa = kappa * (1e-6) * (60*60*24*365);

% Spacing
dz = z(2)-z(1);

% Calculate dt based on spacing for stability
dt = dz^2/(2*kappa);

% Create Initial temperature using as the adiabat
T = Tp + gradT*z;
n = length(T);

% initialize time at t = 0
t = 0;

% Set boundary conditions
T(1) = To;

% Loop to get desired age
Tnew = T;

const= kappa*dt/dz^2;

while (t <=age)
    
    
    %Tnew(2:n-1) = T(2:n-1) + const*(T(2:n)-2*
    for i = 2:(n-1)
        Tnew(i) = T(i) + const*(T(i+1)-2*T(i)+T(i-1));
    end
    
    Tnew(1) = To;
    Tnew(end) = T(end);    
    T = Tnew;
    
    t = t + dt;
end
