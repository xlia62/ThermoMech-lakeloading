function z = marine_deposition(x,z,qs,total_time,Ks,Lcrit)

% z is the bathimetry function
% dx is the node spacing
% qs is the sediment flux from continent at left
% Ks is the sediment diffusion coefficient
% dt is the total diffusion time

delta_t = 100;  % years
% Marine depostion
%qs = 2000;  % m^3
n = length(z);

if n == 1
    return
end

depth = -z;
dx = x(2) - x(1);

% Sediment per node
qs = qs/dx;

% How many nodes in Lcrit, restricted to size of basin
ns = round(Lcrit/dx); % Nodes to distribute to
nn = min(ns,n);
n0 = n-nn;  % Remaining nodes to fill if need be

% Sediment thickness per node
ds = qs/nn;
j = 1:nn;

% Distribute sediment to each node
z(j) = z(j) + ds;

% Check for overfill, and if present, enter loop to redistribute
k = z > 0;
ds = sum(z(k));
z(k) = 0;

while ds > 0 && n0 >= 0
    
    k = z(j) < 0;  % Remaining nodes in Lcrit left to fill
    if any(k)
        nn = sum(k);
        ds = ds/nn;
        z(j(k)) = z(j(k)) + ds;
        k = z > 0;
        if any(k)
            ds = sum(z(k));
            z(k) = 0;
        else
            ds = 0;
        end
    elseif n0==0
        n0 = -1;
    else % Need to add more nodes to list
        nn = min(ns,n0);
        j = j + nn;
        n0 = n0 - nn;
    end
    
end

if ds > 0
    disp('Ocean basin is completely filled')
end

%
% % Sediment per node
% qs = qs/dx;
%
%
%
%
%
%
%
% i = 1;
% ha = depth(i);
% ds = min(qs,ha);
%
% z(i) = z(i) + ds;
% qs = qs - ds;
% while qs > 0 && i < n
%     i = i + 1;
%     ha = depth(i);
%     ds = min(qs,ha);
%     z(i) = z(i) + ds;
%     qs = qs - ds;
% end
%
% if qs > 0
%     % Filled ocean, no need to apply diffusion of sediments
%     disp('Ocean basin completely filled')
% else
if Ks > 0
    % Apply diffusion of sediment with Neumann boundary conditions at source
    % node and a Dirichlet boundary condition at most the critical distance
    k = x <= x(1)+10*Lcrit;
    jz = z(k);
    
    nn = length(jz);
    
    % Apply diffusion for sediment distribution, using implicit scheme
    % solving the system of equations
    
    % Assemble the rhs vector
    rhs = jz(:);
    
    % Assamble sparse coefficient matrix off-diagonals
    alpha = -Ks*delta_t/dx^2;
    A2 = sparse(2:nn,1:nn-1,alpha,nn,nn);
    
    % Assamble sparse coefficient matrix diagonal and sum with off-diagonals
    A1 = sparse(1:nn,1:nn,1-2*alpha,nn,nn);
    A = A2 + A1 + A2';
    
    % Apply the boundary conditions
    % Left BC is zero flux, i.e. jz(1) = jz(2)
    % Right BC is fixed to initial value of jz(end)
    A(1,1) = 1-2*alpha;
    A(1,2) = 2*alpha;
    A(nn,nn) = 1;
    A(nn,nn-1) = 0;
    
    % Solve, however we need to do this in small timesteps because of the
    % high diffusivity we will use the courant condition even though this
    % is an implicit scheme.
    t = 0;
    nsteps = round(total_time/delta_t);  % A slight error in time steps
    for i = 1:nsteps
        rhs = A\rhs;
    end
    jz = rhs';
    z(k) = jz;
    
    % Account for numerical error and bring all values to sea level that
    % are witnin 0.1 m
    k = z > -0.1;
    z(k) = 0;
end