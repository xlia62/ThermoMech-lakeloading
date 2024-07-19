function [jX,jH,jB,jS] = ErosionDeposition1D(delta_t,total_time,jX,jH,jB,Ld,G,Kf,Kd)

lc = 1e-2; % Lacie's coefficient, W = c*Qw^0.5, note that Theunisen and Huismans (2019)
           % lc is 100x that of Yuan et al., (2019) based on Montgomery &
           % Gran (2001). The discrepancy is taken care of by the
           % difference in deposition coefficient G (Yuan et al., 2019) and
           % 1/alpha (Theunisen and Huismans, 2019) for p = 1 m/yr.
h = 5/3;   % Hack's law (1957) original exponent
m = 0.5;   % Stream power law constants
n = 1;     % Stream power law constants <---- Important!!! leave fixed to 1 for now (implied in code)
Ka = 6.7;  % Hack's law constant
p = 1.0;   % Average precipitation, m/yr <---- leave at 1 for now, will need to implement
           % different values later

toleranceH = 1e-5;  % Sensitivity of max/min search, known as Minimum Prominence.

% The minimum prominence seams to not work well for minimums near zero
k = abs(jH) < 1e-8;
jH(k) = 0;

% This function was bench-marked againts Kellin Whipple's 1D stream profile
% evolution code (detach_ftfs.m) for p = 1 (implied) and G = 0, Kd = 0, m = 0.5 and n =1. 
           
%Kf = 2e-5; % Erosion (incision) coefficient
%Kd = 0.1;%1.07; % Transport Diffusion or diffusivity, hillslope processes, m^2/yr

% Characteristic length scale distance from ridge at which both hillslope
% and river incision are equal. Not used, but a useful measure
% Le = (Kd/(Kf*Ka^m))^(1/(m*h+1));

% Number of nodes
nn = length(jX);

% Precipitation at nodes, Not used, just a place holder for precipiation
% variability
%jP = p + zeros(size(jX));  % m/yr

% Sed flux storage
jS = zeros(size(jH));

% For now using a uniform grid.
dx = jX(2)-jX(1);

% Assable part of the uniform diffusion coefficient matrix that will be
% completed and used at the end of each time step.
if Kd > 0
    b = -Kd/dx^2;
    % Assemble part of the sparse coefficient matrix
    A2 = sparse(2:nn,1:nn-1,b,nn,nn);
end

t = 0;
while t < total_time
    t = t+delta_t;
        
    % Reduce delta_t
    if t > total_time
       delta_t = delta_t - (t - total_time);
       t = total_time;
    end
    
    % Perform incision only if Kf > 0
    
    if Kf > 0
        % Get the node numbers corresponding to the local maxima and minima
        jmin = find(islocalmin(jH,'FlatSelection','center','MinProminence',toleranceH));
        jmax = find(islocalmax(jH,'FlatSelection','center','MinProminence',toleranceH));
        
        nmax = length(jmax);
        
        % Local index for each stream segment increases upslope. Each segment is
        % described by
        %
        % S.n = local number of nodes
        % S.j = local nodes to global nodes map
        % S.Qm = Discharge at each node, using Hack's law and applying m exponent
        % S.dx = distance between local nodes
        % S.z = local elevation
        %
        
        % Check that we have at least one maximum or one local minimum
        
        if isempty(jmin) && isempty(jmax)
            
            if abs(min(jH)-max(jH)) < toleranceH
                % The topography is completely flat; no erosion or depostion
                nseg = 0;
                
            elseif jH(1) < jH(end)
                % Uphill to the right, could also be step function
                nseg = 1;
                
                Sright.x = jX(end) - jX;
                Sright.n = length(Sright.x);
                Sright.Qw = Ka*p*Sright.x.^h; 
                Sright.Qw(end) = 0.1; % Assume ridge node has some discharge, privents division by zero.
                
                % Assemble the local segment
                S(1).j = 1:nn; % Local nodes
                S(1).n = nn;
                S(1).r = [1 1:nn-1];    % Base level at left most node
                S(1).Dij = get_donorsmat(S(1).r,S(1).n);
                S(1).stack = 1:nn;      % Stack is simple here
                S(1).z = jH;
                S(1).Qw = Sright.Qw;
                
            else
                % Uphill to the left
                nseg = 1;
                
                Sleft.x = jX - jX(1);  % Ensure ridge is at 0.
                Sleft.n = length(Sleft.x);
                Sleft.Qw = Ka*p*Sleft.x.^h;  % Use Hack's Law
                Sleft.Qw(1) = 0.1;           
                
                % Assemble the local segment
                S(1).j = 1:nn; % Local nodes
                S(1).n = nn;
                S(1).r = [2:nn nn];    % Base level at right most node
                S(1).Dij = get_donorsmat(S(1).r,S(1).n);
                S(1).stack = nn:-1:1;      % Stack is simple here
                S(1).z = jH;
                S(1).Qw = Sleft.Qw;
                
            end
            
        else  % There must be more than one stream segment
            
            % Note the boundary nodes are not caught as local max or min,
            % thus this works.
            
            % How many segments
            if isempty(jmax)
                % There is one local minimum and thus 1 catchment segment
                nseg = 1;
                
                % Left branch (uphill to left)
                Sleft.x = jX(1:jmin) - jX(1);  % Ensure ridge is at 0.
                Sleft.n = length(Sleft.x);
                Sleft.Qw = p*Ka*Sleft.x.^h;  % Use Hack's Law
                Sleft.Qw(1) = 0.1;
                
                % Right branch (uphill to the right)
                Sright.x = jX(end) - jX(jmin:end);
                Sright.n = length(Sright.x);
                Sright.Qw = p*Ka*Sright.x.^h;
                Sright.Qw(end) = 0.1;
                
                S(1).j = 1:nn;
                S(1).n = length(S(1).j);
                S(1).r = [2:jmin jmin jmin:nn-1];
                S(1).Dij = get_donorsmat(S(1).r,S(1).n);
                S(1).stack = [jmin:-1:1 jmin+1:nn];
                
                S(1).z = jH;
                S(1).sed = zeros(1,S(1).n);
                S(1).Qw = [Sleft.Qw(1:end-1) Sleft.Qw(end)+Sright.Qw(1) Sright.Qw(2:end)];                
                
            elseif isempty(jmin)
                % There is one local maximum and thus 2 segments.
                nseg = 2;
                
                % First, uphill to the right
                % The left most node will be base level
                
                Sright.x = jX(jmax) - jX(1:jmax);
                Sright.n = length(Sright.x);
                Sright.Qw = p*Ka*Sright.x.^h;
                Sright.Qw(end) = 0.1;
                
                S(1).j = 1:jmax;
                S(1).n = length(S(1).j);
                S(1).r = [1 1:jmax-1];    % Base level at left node
                S(1).Dij = get_donorsmat(S(1).r,S(1).n);
                S(1).stack = 1:jmax;
                S(1).z = jH(S(1).j);
                S(1).Qw = Sright.Qw;
                
                % Uphill to the left
                % The right most node will be base level
                Sleft.x = jX(jmax:end) - jX(jmax);  % Ensure ridge is at 0.
                Sleft.n = length(Sleft.x);
                Sleft.Qw = p*Ka*Sleft.x.^h;
                Sleft.Qw(1) = 0.1;
                
                S(2).j = jmax:nn;
                S(2).n = length(S(2).j);
                S(2).r = [2:S(2).n S(2).n];  % Base level at right node
                S(2).Dij = get_donorsmat(S(2).r,S(2).n);
                S(2).stack = [nn:-1:jmax] - jmax + 1;  % Local coordinate stack
                S(2).z = jH(jmax:end);
                S(2).Qw = Sleft.Qw;
                
            else
                % There are more than one local minima, each local minimum will
                % form a single catchemnt.
                
                % There is a special case, the boundary nodes are not
                % caught if they are the minimum or maximum
                rmin_left = false;
                if jmax(1) > jmin(1)
                    if jH(1) > jH(jmin(1))
                        jmax  = [ 1 jmax];
                    end
                elseif jmin(1) > jmax(1)
                    if jH(1) < jH(jmax(1)) 
                        jmin = [1 jmin];
                        rmin_left = true;
                    end
                end

                rm_right = false;
                if jmin(end) > jmax(end)
                   if jH(end) > jH(jmin(end))
                       jmax = [jmax nn];
                   end
                elseif jmax(end) > jmin(end)
                    if jH(end) < jH(jmax(end))
                        jmin = [jmin nn];
                        rm_right = true;
                    end
                end
                
                nseg = 2*nmax + 2; %Initial number of segments, could increase by 2
                
                %             % Preallocae space
                %             S = struct('j', cell(1, nseg), 'n', cell(1, nseg),'x', cell(1, nseg),...
                %                 'dx', cell(1, nseg),'z', cell(1, nseg),'p', cell(1, nseg),...
                %                 'KQm', cell(1, nseg));
                %
                % Number of segments
                c = 0;
                
                % First, consider the case that there is a segment to left of
                % the first catchment
                if jmin(1) < jmax(1)
                    c = c + 1;
                    
                    lmax = jmax(1);
                    % Slope up to the right
                    Sright.x = jX(lmax) - jX(1:lmax);
                    Sright.n = length(Sright.x);
                    Sright.Qw = p*Ka*Sright.x.^h;
                    Sright.Qw(end) = 0.1;
                    
                    S(c).j = 1:lmax;
                    S(c).n = length(S(c).j);
                    S(c).r = [1 1:S(c).n-1];    % Base level at left node
                    S(c).Dij = get_donorsmat(S(c).r,S(c).n);
                    S(c).stack = 1:S(c).n;                    
                    S(c).z = jH(S(c).j);
                    S(c).Qw = Sright.Qw;
                    
                    % Remove the local boundary minimum
                    if rmin_left
                        jmin(1) = [];
                    end
                end
                
                % Now consider the case of a segment to right of the last
                % castchment
                if jmax(end) < jmin(end)
                    c = c + 1;
                    
                    lmax = jmax(end);
                    % Slope up to the left
                    
                    Sleft.x = jX(lmax:end)-jX(lmax);  % Ensure ridge is at 0.
                    Sleft.n = length(Sleft.x);
                    Sleft.Qw = p*Ka*Sleft.x.^h;
                    Sleft.Qw(1) = 0.1;
                    
                    S(c).j = lmax:nn;
                    S(c).n = length(S(c).j);
                    S(c).r = [2:S(c).n S(c).n];    % Base level at right node, no ghost base
                    S(c).Dij = get_donorsmat(S(c).r,S(c).n);
                    S(c).stack = S(c).n:-1:1; 
                    S(c).z = jH(S(c).j);
                    S(c).Qw = Sleft.Qw;
                    
                    % Remove the local boundary minimum
                    if rm_right
                        jmin(end) = [];
                    end
                end
                
%                 % Fix up the ends (if flat before descent to minimum)
%                 if jmin(1) < jmax(1) && jH(jmin(1)) < jH(1)
%                     jmax = [1 jmax];
%                 end
%                 
%                 % if flat after last ascent from mimimum
%                 if jmin(end) > jmax(end) && jH(jmin(end)) < jH(end)
%                     jmax(end+1) = nn;
%                 end
                
                for i = 1:length(jmin)
                    c = c + 1;
                    
                    % Global nodes
                    gnodes = jmax(i):jmax(i+1);
                    gnn = length(gnodes);
                    
                    % Left branch (uphill to left)
                    Sleft.x = jX(jmax(i):jmin(i)) - jX(jmax(i));  % Ensure ridge is at 0.
                    Sleft.n = length(Sleft.x);
                    Sleft.Qw = p*Ka*Sleft.x.^h;
                    Sleft.Qw(1) = 0.1;
                    
                    % Right branch (uphill to the right)
                    Sright.x = jX(jmax(i+1)) - jX(jmin(i):jmax(i+1));
                    Sright.n = length(Sright.x);
                    Sright.Qw = p*Ka*Sright.x.^h;
                    Sright.Qw(end) = 0.1;
                  
                    % mimum in local nodes coordinate
                    lmin = jmin(i) - jmax(i) + 1;
                    
                    S(c).j = gnodes;
                    S(c).n = gnn;
                    S(c).r = [2:lmin lmin lmin:S(c).n-1];
                    S(c).Dij = get_donorsmat(S(c).r,S(c).n);
                    S(c).stack = [lmin:-1:1 lmin+1:S(c).n];
                    S(c).z = jH(gnodes);
                    S(c).Qw = [Sleft.Qw(1:end-1) Sleft.Qw(end)+Sright.Qw(1) Sright.Qw(2:end)];
                    
                end
                
                nseg = c;
                
            end
            
        end
        
        % Following Braun and Willet (2013), using implict scheme
        for k = 1:nseg
            
            Zinit = S(k).z;
            Znew = S(k).z;

            % Loop over nodes from just above base level (local min) to just before ridge
            % (local max)  Braun and Willet (2013) EQ 22.
            tol = 0.01;
            if G > 0
                
                err = 1;
                F = Kf*delta_t*S(k).Qw.^m / dx;
                Gfact = G*lc./S(k).Qw;  % note ridges are assigned a Qw that is very small
                maxit = 40;
                numit = 0;
                b = zeros(1,S(k).n);
                ups_nodes = cell(1,S(k).n);
                
                % To speed up the loop, we perform the first itteration
                % outside of the loop to calculate the upstream nodes.

                while err > tol && numit < maxit
                    numit = numit + 1;
                    Zold = Znew;
                    
                    for j = 2:S(k).n
                        i = S(k).stack(j);
                        r = S(k).r(i);
                        if numit == 1
                            ups_nodes{i} = find_upstream(i,S(k).Dij,S(k).n);
                            b(i) = S(k).z(i) + Gfact(i)*sum(S(k).z(ups_nodes{i}).*S(k).Qw(ups_nodes{i}));
                        end
                        Znew(i) = b(i) - Gfact(i)*sum(Zold(ups_nodes{i}).*S(k).Qw(ups_nodes{i})) + F(i)*Znew(r);
                        Znew(i) = Znew(i)/(1+F(i));
                    end
                    err = max(abs(Znew - Zold));
                                     
                end
                if numit == maxit && err > tol
                    warning('max iterations reached')
                    disp(num2str(err))
                end
                
            else
                
                for j = 2:S(k).n
                    i = S(k).stack(j);
                    r = S(k).r(i);
                    F = Kf*delta_t*S(k).Qw(i)^m / dx;
                    Znew(i) = (Znew(i) + Znew(r)*F)/(1 + F);
                end
                % Now map the local segment back into the global topography profile
            end
            
            % Account for numerical round-off errors in flat regions
            j = abs(Zinit-Znew) < tol;
            Znew(j) = Zinit(j);
            
            % Calculate ammount of sediment to deposit
            % This is already an underestimate, because it only takes into account
            % elevation change along the stream profile and not it's tributaries
            
            % Total sediment discharge per unit width at base level node
            % since base-level node is fixed it is quick to calculate
            qs = dx*sum(Zinit - Znew);
            
            % Base node
            ir = S(k).stack(1);
            
            % Accumulate and store the base level sediment flux and map into global array
            jS(S(k).j(ir)) = qs + jS(S(k).j(ir));
            
            % Sediment agrradation at base level node that are interior
            if Ld                
        
                % We want to distribute the sediment over each node,
                % Total sediment pernode
                qs = qs/dx;
                
                % If we are at edges of model space, no sediment is deposited,
                % all flow out of the model
                if ir ~= 1 || ir ~= S(k).n
                    
                    % Calculate the sill height
                    sill = min(Znew(1),Znew(end));
                    % Set the starting sediment increment according to
                    ds = min(qs,sill-Znew(ir))/2; % At most fill half-way
                    Zi = Znew(ir) + ds;  % Horizon level
                    maxit = 100;
                    c = 0;
                    
                    while qs > 0.1 && c < maxit && Zi < sill
                        c = c + 1;
                        dZ = Zi - Znew;
                        i = dZ < 0;  % If negative, horizon is below a physical node
                        dZ(i) = 0;   % set dh to zero at those positions
                        % Subtract deposited sediment from total qs (note, dx is constant)
                        if qs >= sum(dZ)
                            % Increment elevation to account for deposition
                            Znew = Znew + dZ;
                            % Remaining sediment to distribute
                            qs = qs - sum(dZ);
                            % Calculate sediment to redistribute
                            ds = min(qs,sill-Znew(ir))/2; % At most fill half-way
                            Zi = Znew(ir) + ds;  % New horizon level
                        else % Reduce the amount of sediment to redistribute
                            ns = sum(i); % number of nodes to distribute to
                            ds = ds/ns;
                            Zi = Znew(ir) + ds;
                        end
                        
                    end
                end  % end accredation loop
                
            end
            % Put segment back into global array
            jH(S(k).j) = Znew;
            
        end  % end segment
        
    end % end incision
    
    % Apply diffusion for hillslope processes, using implicit scheme
    % solving the system of equations
    
    if Kd > 0
        % Assemble the rhs vector
        rhs = jH'/delta_t;
        
        % Finish constructing the diffusion coefficient matrix taking into
        % account possible change in delta_t
        a = (1/delta_t + 2*Kd/dx^2);
        A1 = sparse(1:nn,1:nn,a,nn,nn);
        A = A2 + A1 + A2';
        
        % Apply the boundary conditions
        A(1,1) = 1/delta_t;
        A(1,2) = 0;
        A(nn,nn) = 1/delta_t;
        A(nn,nn-1) = 0;
        
        % Solve
        jH = A\rhs;
        jH = jH';
    end
    
    % Generate a zero
    k = abs(jH) < 1e-8;
    jH(k) = 0;
    
    % Update the basement
    k = jH - jB < 0;
    jB(k) = jH(k);
    
end  % end time loop

% Now save only the boundary sediment flux
jS = jS([1 end]);

end

%% Suporting internal functions
%

function Dij = get_donorsmat(r,np)

% Assemble a donor matrix 
% following Braun and Willet (2013)

% Get list of donor nodes to each node i
d = zeros(1,np);
for i = 1:np
    d(r(i)) = d(r(i)) + 1;
end

% Calculate delta index from the above list
delta = zeros(1,np+1);
delta(np+1) = np+1;
for i = np:-1:1
    delta(i)=delta(i+1)-d(i);
end

% Generate list of donors
D = zeros(1,np);
w = zeros(1,max(r));
for i = 1:np
    ri = r(i);
    D(delta(ri)+w(ri)) = i;
    w(ri) = w(ri) + 1;
end

% Generate donor information (connectivity)
Dij = zeros(3,np);

for i = 1:np
    for j = 1:delta(i+1)-delta(i)        
        Dij(j,i) = D(delta(i)+j-1);
    end
end

end % end get_donorsmat

function ups_nodes = find_upstream(node,D,nmax)

ups_nodes = NaN(1,nmax);

c = 0;
doner = D(:,node)';
D1 = D(1,:);
for k = 1:length(doner)
    i = doner(k);    
    while i > 0 && i ~= node
        c = c +1;
        ups_nodes(c) = i;
        i = D1(i);
    end
end

ups_nodes = ups_nodes(1:c);

end % end find_upstream

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
