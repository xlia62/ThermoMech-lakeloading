function [gridt,H,E] = surface_topography2(gridx,gridy,xstp1,ystp1,gridcx,gridcy,xstpc1,ystpc1,vx,vy,gridt,H,E,tstp,LEMpar,timestep)
% Surface Topography 2
%
% V1.0 Robert Moucha, Syracuse University
% September 17, 2020
%
% Applies Peice-wise stream power lower with hill-slope processes and
% deposition of sediments


% Set velocity initially to zero
gridt(4:1:6,:)=0;

tnum = size(gridt,2);
xnum = length(gridx);
ynum = length(gridy);

%--------------------------------------------------------------------------
%% Defining advection velocity at topography points
%--------------------------------------------------------------------------
for i=1:1:tnum
    
    % Check topography nodes inside the grid
    if (gridt(1,i)>=gridx(1) && gridt(1,i)<=gridx(xnum) && gridt(2,i)>=gridy(1) && gridt(2,i)<=gridy(ynum))
        %  xn    V(xn,yn)--------------------V(xn+1,yn)
        %           ?           ^                  ?
        %           ?           ?                  ?
        %           ?          dy                  ?
        %           ?           ?                  ?
        %           ?           v                  ?
        %           ?<----dx--->o Mrho(xm,ym)       ?
        %           ?                              ?
        %           ?                              ?
        %  xn+1  V(xn,yn+1)-------------------V(xn+1,yn+1)
        
        % Define indexes for upper left BASIC node in the cell where the topograhy node is
        % using bisection
        xcur=gridt(1,i);
        ycur=gridt(2,i);
        % Find horizontal index
        xnmin=1;
        xnmax=xnum;
        while ((xnmax-xnmin)>1)
            % !!! SUBTRACT 0.5 since int16(0.5)=1
            xn=double(int16((xnmax+xnmin)./2-0.5));
            if(gridx(xn)>xcur)
                xnmax=xn;
            else
                xnmin=xn;
            end
        end
        % Check horizontal index
        if (xnmin<1)
            xnmin=1;
        end
        if (xnmin>xnum-1)
            xnmin=xnum-1;
        end
        
        % Find vertical index
        ynmin=1;
        ynmax=ynum;
        while ((ynmax-ynmin)>1)
            % !!! SUBTRACT 0.5 since int16(0.5)=1
            yn=double(int16((ynmax+ynmin)./2-0.5));
            if(gridy(yn)>ycur)
                ynmax=yn;
            else
                ynmin=yn;
            end
        end
        % Check vertical index
        if (ynmin<1)
            ynmin=1;
        end
        if (ynmin>ynum-1)
            ynmin=ynum-1;
        end
        
        % Define indexes for upper left node in the Vx-cell where topography node is
        % Horizontal Vx index
        xn=xnmin;
        % Vertical Vx index
        yn=ynmin;
        if(ycur>gridcy(yn+1))
            yn=yn+1;
        end
        if (yn>ynum)
            yn=ynum;
        end
        
        % Define and check normalized distances from topography node to the upper left VX-node;
        dx = (xcur-gridx(xn))./xstp1(xn);
        dy = (ycur-gridcy(yn))./ystpc1(yn);
        
        % Calculate topography point velocity from four surrounding Vx nodes
        gridt(4,i)=gridt(4,i)+(1.0-dx).*(1.0-dy).*vx(yn,xn);
        gridt(4,i)=gridt(4,i)+(1.0-dx).*dy.*vx(yn+1,xn);
        gridt(4,i)=gridt(4,i)+dx.*(1.0-dy).*vx(yn,xn+1);
        gridt(4,i)=gridt(4,i)+dx.*dy.*vx(yn+1,xn+1);
        
        % Define indexes for upper left node in the VY-cell where the topography node is
        % Vertical Vy index
        yn=ynmin;
        % Horizontal Vy index
        xn=xnmin;
        if(xcur>gridcx(xn+1))
            xn=xn+1;
        end
        if (xn>xnum)
            xn=xnum;
        end
        
        % Define and check normalized distances from topography node to the upper left VX-node;
        dx=(xcur-gridcx(xn))./xstpc1(xn);
        dy=(ycur-gridy(yn))./ystp1(yn);
        
        % Calculate topography node velocity from four surrounding Vy nodes
        gridt(5,i)=gridt(5,i)+(1.0-dx).*(1.0-dy).*vy(yn,xn);
        gridt(5,i)=gridt(5,i)+(1.0-dx).*dy.*vy(yn+1,xn);
        gridt(5,i)=gridt(5,i)+dx.*(1.0-dy).*vy(yn,xn+1);
        gridt(5,i)=gridt(5,i)+dx.*dy.*vy(yn+1,xn+1);
    end
end


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%% Advect topography and apply surface process model

yr2sec = 365.25*24*3600;
dt_max = yr2sec*LEMpar.dt_max;

% Take into account timesteps lower than desired dt_max
dt_max = min(timestep,dt_max);
dt_extra = rem(timestep,dt_max);

ntimes = floor(timestep/dt_max);

for ijk = 1:ntimes
    
    if ijk == ntimes
        dt_max = dt_max + dt_extra;
    end
    
    %--------------------------------------------------------------------------
    %% Advect topography vertically
    %--------------------------------------------------------------------------
    
    gridt(2,:) = gridt(2,:) + gridt(5,:)*dt_max;
    
    % Advect erosional surface
    E = E + gridt(5,:)*dt_max;
    
    % Advect the horizons vertically
    nH = size(H,1);
    for i = 1:nH
        H(i,:) = H(i,:) + gridt(5,:)*dt_max;
    end
    
    %--------------------------------------------------------------------------
    %% Advect topography horizontally
    %--------------------------------------------------------------------------
    
    % Define maximal horizontal velocity at topography nodes
    vxmax=max(abs(gridt(4,:)));
    
    % Defining topography advection timestep
    dt = dt_max;
    nt = 1;
    if(vxmax>0)
        dt = tstp/vxmax;
        if (dt<dt_max)
            nt = double(int16(dt_max./dt-0.5))+1;
            dt = dt_max/nt;
        else
            dt = dt_max;
        end
    end
    
    % Defining FCT parameter MU
    mu=1/8;
    
    % Advect topography, horizons and erosinal surface with FCT
    for t=1:1:nt
        % Step 0: Set new profile
        gridt(3,:)=gridt(2,:);
        tmpH = H;
        tmpE = E;
        % Step 1: Transport+numerical diffusion stage
        for i=2:1:tnum-1
            % Defining FCT parameters EPS and NU
            eps=gridt(4,i)*dt/tstp;
            nu=1/8+(eps^2)/2;
            % Change topography
            gridt(3,i)=gridt(2,i)-eps/2*(gridt(2,i+1)-gridt(2,i-1))+nu*(gridt(2,i+1)-2*gridt(2,i)+gridt(2,i-1));
            % Change horizons
            tmpH(:,i)=H(:,i)-eps/2*(H(:,i+1)-H(:,i-1))+nu*(H(:,i+1)-2*H(:,i)+H(:,i-1));
            tmpE(i)=E(i)-eps/2*(E(i+1)-E(i-1))+nu*(E(i+1)-2*E(i)+E(i-1));
        end
        % Step 2: Antidiffusion stage
        % Antidiffusion flow for the first cell
        gridt(6,1)=0;
        tmpH2 = zeros(size(H));
        tmpE2 = zeros(size(E));
        for i=2:1:tnum-2
            % Corrected antidiffusion flow for current cell
            delt0=gridt(3,i)-gridt(3,i-1);
            delt1=gridt(3,i+1)-gridt(3,i);
            delt2=gridt(3,i+2)-gridt(3,i+1);
            s=sign(delt1);
            gridt(6,i)=s*max(0,min(min(s*delt2,s*delt0),mu*abs(delt1)));
            gridt(2,i)=gridt(3,i)-gridt(6,i)+gridt(6,i-1);
            % Correct antidiffusion flow for current horizon cell
            delt0=tmpH(:,i)-tmpH(:,i-1);
            delt1=tmpH(:,i+1)-tmpH(:,i);
            delt2=tmpH(:,i+2)-tmpH(:,i+1);
            s = sign(delt1);
            for k = 1:nH
                tmpH2(k,i)=s(k)*max(0,min(min(s(k)*delt2(k),s(k)*delt0(k)),mu*abs(delt1(k))));
                H(k,i)=tmpH(k,i)-tmpH2(k,i)+tmpH2(k,i-1);
            end
            % Correct antidiffusion flwo for current erosional surface
            delt0=tmpE(i)-tmpE(i-1);
            delt1=tmpE(i+1)-tmpE(i);
            delt2=tmpE(i+2)-tmpE(i+1);
            s = sign(delt1);
            tmpE2(i)=s*max(0,min(min(s*delt2,s*delt0),mu*abs(delt1)));
            E(i)=tmpE(i)-tmpE2(i)+tmpE2(i-1);
        end
    end
    
    %--------------------------------------------------------------------------
    %% Apply the surface processes model
    %--------------------------------------------------------------------------
    
    if LEMpar.apply_surfacepro
        
        % Run the model only if maximum relief exceeds 20 m
        if max(gridt(2,:)) - min(gridt(2,:)) >= 20
            
            % Since topography is in the model y-space, i.e. basins have a greater
            % value that mountains, we subtract this value and thus set the lowest
            % value in topogrpahy to 0 m elevation. This also accounts for varying
            % sticky air thickness
            offset = max(gridt(2,:));
            % Note that the first horizon is the basement.
            [~,temp_topo,temp_H] = ErosionDeposition1D(dt_max/yr2sec,dt_max/yr2sec,gridt(1,:),offset-gridt(2,:),offset-H(1,:),LEMpar.Ld,LEMpar.Fd,LEMpar.Kf,LEMpar.Kd);
            
            % Bring the topography back into the model y-space
            gridt(2,:) = offset-temp_topo;
            H(1,:) = offset-temp_H;
        end
        
        % Make sure that the horizons do not cross the topography and reflect the
        % basement
        for i = 1:nH
            k = gridt(2,:) > H(i,:);
            H(i,k) = gridt(2,k);
        end
        
    end
    
end

end %function

