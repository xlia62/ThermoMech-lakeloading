% Version 3.1 August 15, 2017
% Rob Moucha, Syracuse University
%
% Eclogite Phase Change -- Equilibrium and Kinetic based with density
%                          change
%                          Fixed equation of state for phase change
%
% Mantle Transizition Phase Changes -- Equilibrium based with density 
%                                      change and latent heat 
%
% Version 3.0 March 14, 2017
% Rob Moucha, Syracuse University
%
% !!!! Changed timing of ouptut from iterations steps to time steps !!!
% - Initial Eclogite Phase Transition
%
% Version 2.1 Feb 7, 2017
% Rob Moucha, Syracuse University
% Siobhan Campbell, Syracuse University
% Prasanna Gunawardana, Syracuse University
%
% - Bug fixes, melt and fluid weakening corrected
%
% Version 2.0 October 1, 2016
% Rob Moucha, Syracuse University
% Siobhan Campbell, Syracuse University
%
% - New marker boundary condition that recycles markers from one side of
%   model as they exit to the other side
%
% - Bug fixes, restart enabled
% 
% Version 1.0 Release 5.0  - July 5, 2016
% Rob Moucha, Syracuse University
%
% - Include melt generation and extraction (not availble in release 4)
%

% Initialize the computational grid variables
% Density, viscosity, shear modulus, temperature, thermal conductivity, RHO*Cp arrays
       
if ~restart
    G.etas1 = zeros(ynum,xnum);       % Viscosity for shear stress
    G.etan1 = zeros(ynum-1,xnum-1);   % Viscosity for normal stress
    G.mus1 = zeros(ynum,xnum);        % Shear modulus for shear stress
    G.mun1 = zeros(ynum-1,xnum-1);    % Shear modulus for normal stress
    G.sxy1 = zeros(ynum,xnum);        % Shear stress
    G.sxx1 = zeros(ynum-1,xnum-1);    % Normal stress
    G.rho1 = zeros(ynum,xnum);        % Density
    G.tk1 = zeros(ynum,xnum);         % Old temperature
    G.tk2 = zeros(ynum,xnum);                        % New temperature
    G.rhocp1 = zeros(ynum,xnum);      % RHO*Cp (for temperature equation)
    G.kt1 = zeros(ynum,xnum);         % Thermal conductivity
    G.hr1 = zeros(ynum,xnum);         % Radiogenic heat production
    G.ha1 = zeros(ynum,xnum);         % Adiabatic heat production/consuming
    G.hph = zeros(ynum,xnum);         % Latent heat from phase changes
    G.he = zeros(ynum,xnum);          % Adiabat/Shear heating efficiency
    G.vx1 = zeros(ynum+1,xnum);
    G.vy1 = zeros(ynum,xnum+1);
    G.pr1 = zeros(ynum-1,xnum-1);
end

G.wtnodes=zeros(ynum,xnum);       % weigths for basic nodes
G.wtetas=zeros(ynum,xnum);        % weights for etas
G.wtetan=zeros(ynum-1,xnum-1);    % weights for etan

TXN = zeros(1,marknum);           % topography marker index

if moving_grid
    gridx0 = G.gridx(1);          % used in moving_grid and topography
    gridxN = G.gridx(end);
end

if movebound  % keep track of internal boundary position
    ix = G.gridx(bintern(1));
end

if ~restart
    % (Re-) Initialize time, s
    timesum = 0;
    startstep = 1;
end

markers_added = 0;
markers_removed = 0;

% Keeping track of when to inject markers and boundaries
bxl = 0;
bxr = 0;
cx = 0;

timestep = Par.timestep/3;

%% Main Time cycle
for ntimestep = startstep:Par.stepmax    
    
    if timesum/yr2sec > 1e6
        disp('here')
    end
    comp_time(1) = 0;  tic
    disp('---------------------------------')
    disp(['Processing time step: ',num2str(ntimestep),' of ',num2str(Par.stepmax)])
    
    % Defining viscoelastic timestep
%     if ntimestep < 3 || restart 
%         timestep = 500*yr2sec;  % First two time steps are 500 years, helps equilibration
%     else
        timestep = min(timestep*3,Par.timemax); % viscoelastic timestep
%     end
%     
    % Plastic yeilding mark
    plastyn=0;
    
    % New markers
    newmarkers = 0;
    
    % Additional markers do to melt extraction
    mm1plus = 0;   
    
    if ntimestep == 1 || moving_grid || restart
        % Computing grid steps for basic nodes
        xstp1 = gridx(2:xnum) - gridx(1:xnum-1);
        ystp1 = gridy(2:ynum) - gridy(1:ynum-1);
        
        % Computing grids and grid steps for Vx, Vy nodes
        % Horizontal (for Vy)
        gridcx=zeros(xnum+1,1);
        xstpc1=zeros(xnum,1);
        % Vertical (for Vx)
        gridcy=zeros(ynum+1,1);
        ystpc1=zeros(ynum,1);
        % First and last nodes and steps (for external nodes)
        % Horizontal (for Vy)
        gridcx(1)=gridx(1)-xstp1(1)/2;
        xstpc1(1)=xstp1(1);
        gridcx(xnum+1)=gridx(xnum)+xstp1(xnum-1)/2;
        xstpc1(xnum)=xstp1(xnum-1);
        % Vertical (for Vx)
        gridcy(1)=gridy(1)-ystp1(1)/2;
        ystpc1(1)=ystp1(1);
        gridcy(ynum+1)=gridy(ynum)+ystp1(ynum-1)/2;
        ystpc1(ynum)=ystp1(ynum-1);
        
        % Internal nodes
        gridcx(2:xnum) = (gridx(2:xnum) + gridx(1:xnum-1))/2;
        gridcy(2:ynum) = (gridy(2:ynum) + gridy(1:ynum-1))/2;
        
        % Internal grid steps
        xstpc1(2:xnum-1) = (gridx(3:xnum) - gridx(1:xnum-2))/2;
        ystpc1(2:ynum-1) = (gridy(3:ynum) - gridy(1:ynum-2))/2;
    end
    
    % Check and delete markers outside the grid
    k = MX < gridx(1) | MX > gridx(xnum) | MY < gridy(1) | MY > gridy(ynum);
    markers_removed = markers_removed + sum(k);
    MX(k) = [];
    MY(k) = [];
    MXN(k) = [];
    MYN(k) = [];
    MCXN(k) = [];
    MCYN(k) = [];
    MI(k) = [];
    MTK(k) = [];
    MRAT(k) = [];
    MGII(k) = [];
    MBII(k) = [];
    MPR(k) = [];
    MSXX(k) = [];
    MSXY(k) = [];
    META(k) = [];
    MEXX(k) = [];
    MEXY(k) = [];
    MXM(k) = [];
    MEXT(k) = [];
    MEXTC(k) = [];
    TXN(k) = [];
    MLAMBDA(k) = [];
    if RxType
        MECL(k,:) = [];
        MPH410(k,:) = [];
        MPH660(k,:) = [];
    end

    plastic_yield(k) = [];

    % Update the markers
    marknum = length(MI);
    
    % Check for movement of markers across fixed boundaries
    % Make sure that water and air markers don't cross-over into the rocks
    
    if sticky_layer > 0
        
        for mm1 = 1:marknum
            
            % Erosion-sedimentation
            % Find topography node for the marker
            xn=double(int16((MX(mm1)-gridt(1,1))/tstp-0.5))+1;
            if (xn<1)
                xn=1;
            end
            if (xn>tnum-1)
                xn=tnum-1;
            end
            % Save horizontal index
            TXN(mm1)=xn;
            
        end
        
        % Compute relative distance to topography node
        dx = (MX-gridt(1,TXN)')/tstp;
        % Compute topograhy elevation above the marker
        dy = gridt(2,TXN)'.*(1-dx)+gridt(2,TXN+1)'.*dx;
        
        % water/air to sediments transformation
        k = (MI == 1 | MI == 2) & MY > dy;
        MI(k) = 3;   % Change marker type to sediment
        MRAT(k) = 1; % Reset strain rate ratio
        MGII(k) = 0; % Reset strain
        
        % Rocks to air transformation
        k = MI > 1 & MI~=2 & MY < dy;
       
        MI(k) = 1;   % Change marker type
        MRAT(k) = 1; % Reset strain rate ratio
        MGII(k) = 0; % Reset strain
        
        % Air to water transformation and vice versa (can't uplift water)
        if ~isempty(water_lev)
            k = MI == 1 & MY > water_lev;
            MI(k) = 2;
            MRAT(k) = 1; % Reset strain rate ratio
            MGII(k) = 0; % Reset strain
            k = MI == 2 & MY < water_lev;
            MI(k) = 1;
            MRAT(k) = 1; % Reset strain rate ratio
            MGII(k) = 0; % Reset strain
        end
        
%          % Reset temperature to ttop for all sticky_markers
%         k = MI == 2 | MI == 1;
%         MTK(k) = ttop;
        
    end
    
    % Assemble a relative location index for markers inside the grid
    
    if ntimestep == 1 || markmove == 0 || restart || moving_grid
        for mm1 = 1:marknum
            
            %  xn    rho(xn,yn)--------------------rho(xn+1,yn)
            %           ?           ^                  ?
            %           ?           ?                  ?
            %           ?          dy                  ?
            %           ?           ?                  ?
            %           ?           v                  ?
            %           ?<----dx--->o Mrho(xm,ym)       ?
            %           ?                              ?
            %           ?                              ?
            %  xn+1  rho(xn,yn+1)-------------------rho(xn+1,yn+1)
            %
            % Define indexes for upper left node in the cell where the marker is
            % using bisection
            % Find horizontal index
            xnmin=1;
            xnmax=xnum;
            while ((xnmax-xnmin)>1)
                % !!! SUBTRACT 0.5 since int16(0.5)=1
                xn = round((xnmax+xnmin)/2-0.5);
                %xn=double(int16((xnmax+xnmin)./2-0.5));
                if(gridx(xn)>MX(mm1))
                    xnmax=xn;
                else
                    xnmin=xn;
                end
            end
            xn=xnmin;
            % Check horizontal index
            if (xn<1)
                xn=1;
            end
            if (xn>xnum-1)
                xn=xnum-1;
            end
            % Save horizontal index
            MXN(mm1)=xn;
            
            % Now define the index for the central nodes (Presure, SXX, etc..)
            if (MX(mm1) < gridcx(xn+1))
                xn = xn-1;
            end
            if (xn<1)
                xn=1;
            end
            if (xn>xnum-2)
                xn=xnum-2;
            end
            MCXN(mm1) = xn;
            
            % Find vertical index
            ynmin=1;
            ynmax=ynum;
            while ((ynmax-ynmin)>1)
                % !!! SUBTRACT 0.5 since int16(0.5)=1
                yn = round((ynmax+ynmin)./2-0.5);                
                %yn=double(int16((ynmax+ynmin)./2-0.5));
                if(gridy(yn)>MY(mm1))
                    ynmax=yn;
                else
                    ynmin=yn;
                end
            end
            yn=ynmin;
            % Check vertical index
            if (yn<1)
                yn=1;
            end
            if (yn>ynum-1)
                yn=ynum-1;
            end
            % Save Vertical index
            MYN(mm1)=yn;
            
            % Now define the index for the central nodes (Presure, SXX, etc..)
            if (MY(mm1) < gridcy(yn+1))
                yn = yn-1;
            end
            if(yn<1)
                yn = 1;
            end
            if(yn>ynum-2)
                yn = ynum-2;
            end
            MCYN(mm1) = yn;
            
        end
    end
   
    comp_time(2) = toc;
    
    % Interpolating parameters from markers to nodes on computational grid
    
    % Define normalized distances from marker to the upper left node;
    
    MDX = (MX - gridx(MXN))./xstp1(MXN);
    MDY = (MY - gridy(MYN))./ystp1(MYN);
    MCDX = (MX - gridcx(MCXN+1))./xstpc1(MCXN+1);
    MCDY = (MY - gridcy(MCYN+1))./ystpc1(MCYN+1);
    
    % Compute density from marker temperature  
    % Account for phase changes if need to
    
    %MRHOCUR = (MRHO(MI,1) + MECL(:,1)*MPHASE(5,1) + MPH410(:,1)*MPHASE(5,3) + MPH660(:,1)*MPHASE(5,4)) ...
    %    .*(1-MRHO(MI,2).*(MTK-273)).*(1+(1-MPH660(:,1).*(1-exp((-MY)/4000e3))).*MRHO(MI,3).*(MPR-1e+5));

    MRHOCUR = MRHO(MI,1).*(1-MRHO(MI,2).*(MTK-273)).*(1+MRHO(MI,3).*(MPR-1e+5));
    
    if RxType
        MRHOCUR = MRHOCUR + MECL(:,1)*MPHASE(5,1) + MPH410(:,1)*MPHASE(5,3) + MPH660(:,1)*MPHASE(5,4);
    end

    % Account for compositional density, i.e. depletion
    MRHOCUR = MRHOCUR + MRHO(MI,5);
    
    % Compute rho*Cp for marker
    MRHOCPCUR = MRHOCUR.*MCP(MI);
    
    % Compute thermal conductivity from marker temperature
    % Rock thermal conductivity (W/m/K): k=k0+a/(T+77)
    MKTCUR = (MKT(MI,1) + MKT(MI,2)./(MTK+77)).*(1 + MKT(MI,3).*MPR);
    
    % Compute adiabatic heating term (alp*T) will be multiplied in heat equation (DP/Dt)
    MHACUR = MRHO(MI,2).*MTK;
    
    comp_time(3) = toc;
    
    % Melting of rocks,
    if (melting && ntimestep>1) % For best results start at 2, for comparison at 0
        
        % Reset melt extract amount
        MEXT(:) = 0;
        
        % Reset fluid weakening
        MLAMBDA(:) = 1;
        
        % Caclulate melt fraction
        [XMELT,~] = Melt_fraction(MPR,MTK,MI);
        
        if extract_melt
            % Total amount of melt
            % Keep track of cummulative melt
            % You can only extract at most 100% melt, so we must account for how
            % much we have already extracted.
            MXM = XMELT - MEXTC;
            
            % Loop over markers to extract melt and form new basalt crust on
            % surface
            k1 = MXM > 0;
            
            % Rock is refractory, means we may have extracted 100% of melt, so
            % it can't contain any melt. If this is negative it also means that
            % we have depeleted mantle
            
            % MXM(~k1) = 0;
            
            % Extract Melt and Create a marker at the surface directly above
            % the melting marker
            
            k2 = MXM > meltmax; % Check for threshold
            MEXT(k2) = MXM(k2) - meltmin; % Amount of extracted melt at this time step
            MEXTC(k2) = MEXTC(k2) + MEXT(k2); % Cummulative extraction of melt
            
            % To save space, only save markers above 200 km depth
            if sum(k2) > 0
                mxmin = min(MX(k2)); mxmax = max(MX(k2));
                mymin = min(MY(k2)); mymax = max(MY(k2));
                
                kout = MX >= mxmin & MX <= mxmax & MY >= mymin & MY <= mymax;
                mxout = MX(kout);
                myout = MY(kout);
                miout = MI(kout);
                mextout = MEXT(kout);
                
                save([outpath,'/meltext_',num2str(ntimestep),'.mat'],'mxout','myout','miout','mextout','timesum')
            end
            
            % Create new markers corresponding to melt extracted directly above
            % at surface, i.e. X-coordinate being the same from which it was
            % extracted, only if there is a sticky_layer
            mm1plus = sum(k2);
            
            if sticky_layer > 0 && mm1plus > 0
                
                marknumnew = length(MI) + (1:mm1plus);
                MX(marknumnew) = MX(k2);
                
                % Use equation of corresponding topography segment
                % to get coordinate of topography on surface
                txn = TXN(k2);
                toposlope = (gridt(2,txn)-gridt(2,txn+1))./(gridt(1,txn)-gridt(1,txn+1));
                MY(marknumnew) = toposlope.*MX(k2)' + (gridt(2,txn) - toposlope.*gridt(1,txn));
                MY(marknumnew) = MY(marknumnew) + 10; % Move slightly below (meters)
                
                TXN(marknumnew) = txn;
                
                % Set the temperature
                MTK(marknumnew) = ttop;
                
                % create basaltic crust by extracting asthenospheric melt
                % Felsic crust can be implemented by looking at the source melt
                MI(marknumnew) = 4;
                
                % Create new marker location indeces, could be spead up if need
                % be
                for i = marknumnew
                    [MXN(i),MYN(i),MCXN(i),MCYN(i)] = locatemarkeringrid(MX(i),MY(i),gridx,gridy,gridcx,gridcy,xnum,ynum);
                end
                
                % Define normalized distances from marker to the upper left nodes
                MDX(marknumnew) = (MX(marknumnew) - gridx(MXN(marknumnew)))./xstp1(MXN(marknumnew));
                MDY(marknumnew) = (MY(marknumnew) - gridy(MYN(marknumnew)))./ystp1(MYN(marknumnew));
                MCDX(marknumnew) = (MX(marknumnew) - gridcx(MCXN(marknumnew)+1))./xstpc1(MCXN(marknumnew)+1);
                MCDY(marknumnew) = (MY(marknumnew) - gridcy(MCYN(marknumnew)+1))./ystpc1(MCYN(marknumnew)+1);
                
                MSXX(marknumnew) = 0;
                MSXY(marknumnew) = 0;
                META(marknumnew) = etamax;
                MEXX(marknumnew) = 0;
                MEXY(marknumnew) = 0;
                MPR(marknumnew) = 1e5;
                MGII(marknumnew) = 0;
                MBII(marknumnew) = 0;
                MRAT(marknumnew) = 1;
                MXM(marknumnew) = 0;
                MEXT(marknumnew) = 0;
                MEXTC(marknumnew) = 0;
                MLAMBDA(marknumnew) = 1;
                if RxType
                    MECL(marknumnew,:) = 0;
                    MECL(marknumnew,2) = 1e4;
                    MPH410(marknumnew,:) = 0;
                    MPH660(marknumnew,:) = 0;
                end
                plastic_yield(marknumnew) = 0;
                
                % Equation of state
                MRHOCUR(marknumnew) = MRHO(MI(marknumnew),1).*(1-MRHO(MI(marknumnew),2).*(MTK(marknumnew)-273)).*(1+MRHO(MI(marknumnew),3).*(MPR(marknumnew)-1e+5));
                % Account for compositional density, i.e. depletion
                MRHOCUR(marknumnew) = MRHOCUR(marknumnew) + MRHO(MI(marknumnew),5);
                % Compute rho*Cp for marker
                MRHOCPCUR(marknumnew) = MRHOCUR(marknumnew).*MCP(MI(marknumnew));
                % Compute thermal conductivity from marker temperature
                % Rock thermal conductivity (W/m/K): k=k0+a/(T+77)
                MKTCUR(marknumnew) = MKT(MI(marknumnew),1) + MKT(MI(marknumnew),2)./(MTK(marknumnew)+77);
                % Compute adiabatic heating term (alp*T) will be multiplied in heat equation (DP/Dt)
                MHACUR(marknumnew) = MRHO(MI(marknumnew),2).*MTK(marknumnew);
                
                % Compute 1/MU values (MU is shear modulus)
                MMUCUR(marknumnew) = 1./MMU(MI(marknumnew));
                
            end
            
            % Locate all markers within a column of 2*dike_dx meters directly
            % above extracted marker and account for the weakening
            % effect of ascending melts (e.g. Vogt et al., 2012)
            if lambda_melt > 0
                k3 = find(k2);  % Focus only rocks that had melt extracted at this time step
                if ~isempty(k3)
                    for i = k3'
                        k = MY < MY(i) & MX < (MX(i)+dike_dx) & MX > (MX(i)-dike_dx) & MI(i) > 2;
                        MLAMBDA(k) = lambda_melt;
                    end
                end
            end
        else
            
            MXM = XMELT;
            k1 = MXM > 0;
        
        end
        
        % Adjust the melt bearing rock paramters
        if any(k1)
                        
            % Reset creep parameters for molten rocks (??? is this
            % warrented ???)
            MRAT(k1) = 1; % Reset strain rate ratio
            MGII(k1) = 0; % Reset plastic strain
            MBII(k1) = 0; % Reset bulk strain
            MLAMBDA(k1) = lambda_melt;  % Weaken rocks containing melt.

            xmelt = XMELT(k1);
            
            % Melt Density
            MRHOCUR(k1) = MRHOCUR(k1).*((1-xmelt)+MRHO(MI(k1),4)./MRHO(MI(k1),1).*xmelt);
            
            % RHO*CP
            MRHOCPCUR(k1) = MRHOCPCUR(k1).*((1-xmelt)+MRHO(MI(mm1),4)./MRHO(MI(k1),1).*xmelt);
            
            % Compute adiabatic heating term (alp*T*DP/Dt)
            MHACUR(k1) = MHACUR(k1).*((1-xmelt)+MRHO(MI(k1),4)./MRHO(MI(k1),1).*xmelt);
            
            % Latent heating: effective adiabatic term, RHOCP
            
            % Melting adiabatic term: alpham=-rho*(dHlat/dP)/T
            % Numerical differentiation
            dp = 1000; % Pressure increment, Pa
            [~, hlat0] = Melt_fraction(MPR(k1)-dp,MTK(k1),MI(k1));
            [~, hlat1] = Melt_fraction(MPR(k1)+dp,MTK(k1),MI(k1));
            MHACUR(k1) = MHACUR(k1) - (hlat1-hlat0)/(2.0*dp);
            
            % Melting heat capacity term: cpm=dHlat/dT
            % Numerical differentiation
            dt = 1.0; % Temperature increment, K
            [~, hlat0] = Melt_fraction(MPR(k1),MTK(k1)-dt,MI(k1));
            [~, hlat1] = Melt_fraction(MPR(k1),MTK(k1)+dt,MI(k1));
            MRHOCPCUR(k1) = MRHOCPCUR(k1) + MRHOCUR(k1).*(hlat1-hlat0)/(2.0*dt);
        end
        
    end

    % Compute the marker rheology
    %[META,MMUCUR] = marker_rheology(MI,MSXX,MSXY,MFLOW,MTK,MEXX,MEXY,MRAT,MPL,MGII,MMU,MPR,ttop,stressmin,etamin,etamax,RGAS,timestep,ntimestep);
    
    % Compute yield stress for each marker
    
    % Plastic strain weakening (Drucker Pager)
    if ntimestep > 1  
        
        % Checking yeilding criterion for strain weakening/hardening
        % C=C0, FI=FI0 for strain<=GAM0
        % C=C0+(C1-C0)/(GAM1-GAM0)*(strain-GAM0), FI=FI0+(FI1-FI0)/(GAM1-GAM0)*(strain-GAM0) for GAM0<strain<GAM1
        % C=C1, FI=FI1 for strain>=GAM1
        
        MCOHES = MPL(MI,1);  % C0
        MFRICT = MPL(MI,3);  % FI0

        %------------------------------------------------------------------
        % This block could be recodede to save memory and flops if needed,
        % Very small hit (0.06 seconds) for 2 million markers
        
        C1 = MPL(MI,2);
        FI1 = MPL(MI,4);

        GAM0 = MPL(MI,5);
        GAM1 = MPL(MI,6);
                
        mid_COHES = MCOHES + (MGII-GAM0).*(C1 - MCOHES)./(GAM1 - GAM0);
        mid_FRICT = MFRICT + (MGII-GAM0).*(FI1 - MFRICT)./(GAM1 - GAM0);
        %------------------------------------------------------------------
        
        % Check for accumalted strain amount
        % Note, here we use the accumalated plastic strain
        k = MGII >= GAM0;
        
        MCOHES(k) = mid_COHES(k);
        MFRICT(k) = mid_FRICT(k);
        
        k = MGII > GAM1;
        
        MCOHES(k) = C1(k);
        MFRICT(k) = FI1(k);
        
        if melting
            siiyield = MCOHES + MFRICT.*MPR.*MLAMBDA;
        else
            siiyield = MCOHES + MFRICT.*MPR;
        end
        
        % Limit the stress to maximum
        siiyield = min(siiyield,MPL(MI,7));
        
    else 
        
        siiyield = MPL(MI,7)*1000;
        
    end
    
    if viscoelastic
        eMMU = MMU;
    else
        eMMU = 1e100*MMU;
    end
    
    % Reset plastic_yield marker
    plastic_yield(:) = -1.0;
    
    
    %checkinputmexrheology(MSXX,MSXY,MFLOW(MI,1:5),MTK,MEXX,MEXY,MRAT,eMMU(MI),MPR,siiyield,plastic_yield,ttop,stressmin,etamin,etamax,RGAS,timestep,ntimestep)
    [META,MMUCUR,MSXX,MSXY,plastic_yield] = mexmarker_rheology(MSXX,MSXY,MFLOW(MI,1:5),MTK,MEXX,MEXY,MRAT,eMMU(MI),MPR,...
        siiyield,plastic_yield,ttop,stressmin,etamin,etamax,RGAS,timestep,ntimestep);    
%    [META,MMUCUR,~,~,plastic_yield] = mexmarker_rheology(MSXX,MSXY,MFLOW(MI,1:5),MTK,MEXX,MEXY,MRAT,eMMU(MI),MPR,...
%        siiyield,plastic_yield,ttop,stressmin,etamin,etamax,RGAS,timestep,ntimestep);    
 
% No idea why plastic strain was reset, should decay with time if no longer
% yielding.

%if ntimestep > 1
    % Plastic strain will accumlate if plastic yield has occured
%    i = plastic_yield < 0;
    %MGII(i) = abs(MGII(i));
    % Else reset plastic strain if no yielding
    %MGII(~i) = -1e-20;
%end
    
    % [META,MMUCUR,MSXX,MSXY] = marker_rheology(MI,MSXX,MSXY,MFLOW,MTK,MEXX,MEXY,MRAT,MPL,MGII,MMU,MPR,ttop,stressmin,etamin,etamax,RGAS,timestep,ntimestep);
    
    % [MSXX,MSXY,META,MMUCUR,plastyn] = marker_rheology_vp(MSXX,MSXY,MEXX,MEXY,MTK,MPR,MRAT,MFLOW(MI,:),MMU(MI),siiyield,ttop,stressmin,etamin,etamax,RGAS,plastyn,ntimestep);
    
    % Because we don't differentiate from sticky_layer in constraining the
    % minimumt viscosity, we need to set the sticky_layer viscosity
    % seperately here
    
    if sticky_layer
        i = MI == 1;
        META(i) = MFLOW(1,2);
        i = MI == 2;
        META(i) = MFLOW(2,2);
    end
    
    %toc
    % Viscosity of partially molten rocks
    % (to be updated with specific rheological rule)
    % After Schmeling (2010)
    if (melting && ntimestep>1) % Only for MXM > 0
        MGII(k1) = 0;  % Make sure no plastic strain in markers with melt
        if fixed_etamelt
            k1 = MXM > 0.1;
            META(k1) = etamelt;
        else
            META(k1) = META(k1).*exp(-28*MXM(k1));
        end
    end
    
    % Reduce viscosity due to melt fraction, but limit to etamelt
    %    META = META.*exp(-28*MXM);
    %    META = max(META,etamelt);

    % Viscosity increase at 660 km
    if RxType
        META = META.*(MPH660(:,1)*(MPHASE(8,4)-1) + 1);
    end

    % Now let's increase the number of markers to account for the new
    % markers. Note that these marker properties will not be used until the
    % next time-step.
    marknum = marknum + mm1plus;
    
    % Computing  Viscosity, density, rock type for nodal points
    % Loop over grid and add properties to 4 surrounding nodes
    
    % Old slower mex routine;
    % Gnew = mexMarker2gridB2(MXN,MYN,MDX,MDY,MRHOCUR,MTK,MKTCUR,MRHOCPCUR,MHR(MI),MHE(MI),MHACUR,META,MMUCUR,MSXX,MSXY,Gnew);
    
    % New mex function
    % Add latent heat production from mantle phase changes to radiogenic
    % heat
    hr = MHR(MI);
    if RxType
        hr = hr + MPH410(:,2) + MPH660(:,2);
    end

    grid_type = G.grid_type;
    vx1old = G.vx1;
    vy1old = G.vy1;
    pr1old = G.pr1;
    
    G = mexMarkers2grid(MXN,MYN,MDX,MDY,MRHOCUR,MTK,MKTCUR,MRHOCPCUR,hr,MHE(MI),MHACUR,META,MMUCUR,MSXX,MSXY,G);
    
    G.grid_type = grid_type;
    
    % Non-mex routine (10 times slower)
    % G = markers2grid(MXN,MYN,MDX,MDY,MRHOCUR,MTK,MKTCUR,MRHOCPCUR,MHR(MI),MHE(MI),MHACUR,META,MMUCUR,MSXX,MSXY,G);
    
    % Flatten density distribution for "air/water" boundary
    % in order to avoid perturbations in the weak layer
    if ~isempty(water_lev)
        for i=1:G.ynum
            for j = 1:G.xnum
                if (G.rho1(i,j)<=1000 && G.gridy(i)<water_lev && G.he(i,j)==0)
                    G.rho1(i,j)=1;
                end
                if (G.rho1(i,j)<1000 && G.gridy(i)>=water_lev && G.he(i,j)==0)
                    G.rho1(i,j)=1000;
                end
            end
        end
    end
    
    comp_time(4) = toc;
    
    % Applying thermal boundary conditions for interpolated temperature
    % Upper, Lower boundaries
    % Upper boundary
    
    G.tk1(1,2:xnum-1) = btopt(2:xnum-1,1)' + btopt(2:xnum-1,2)'.*G.tk1(2,2:xnum-1);
    
    % Lower boundary
    G.tk1(ynum,2:xnum-1) = bbottomt(2:xnum-1,1)' + bbottomt(2:xnum-1,2)'.*G.tk1(ynum-1,2:xnum-1);
    
    % Left boundary
    G.tk1(:,1) = bleftt(:,1) + bleftt(:,2).*G.tk1(:,2);
    % Right boundary
    G.tk1(:,xnum) = brightt(:,1) + brightt(:,2).*G.tk1(:,xnum-1);
    
    % Check for zero viscosity -- indicates not enough markers in cell
    if (min(G.etas1(:))==0)
        disp('Zero viscosity found, interpolating with nearest non-zero')
        k = G.etas1 == 0;
        eta = G.etas1;
        eta(k) = NaN;
        eta = fillmissing(eta,'nearest');
        G.etas1 = eta;
    end
    
    if (min(G.etan1(:))==0)
        disp('Zero viscosity found, interpolating with nearest non-zero')
        k = G.etan1 == 0;
        eta = G.etan1;
        eta(k) = NaN;
        eta = fillmissing(eta,'nearest');
        G.etan1 = eta;
    end
    
    
    % Interpolating initial nodal temperatures back to markers
    % to avoid initial discrepancies between markers and nodes -- done only
    % at the first time_step
 
    if (timesum==0)
        
        % Marker cycle
        for mm1 = 1:marknum
            
            % Markers have not moved, so no need to locate them again
            xn = MXN(mm1);
            yn = MYN(mm1);
            dx = MDX(mm1);
            dy = MDY(mm1);
            
            tkm = (1.0-dx).*(1.0-dy).*G.tk1(yn,xn);
            tkm = tkm + (1.0-dx).*dy.*G.tk1(yn+1,xn);
            tkm = tkm + dx.*(1.0-dy).*G.tk1(yn,xn+1);
            tkm = tkm + dx.*dy.*G.tk1(yn+1,xn+1);
            % Reset marker temperature
            MTK(mm1)=tkm;
            
        end
    end
    
    % Computing viscoelastic (numerical) viscosity and stress
    if viscoelastic
        % Shear stress
        %Viscoelasticity factor
        xelvis = G.etas1./(G.etas1 + timestep*G.mus1);
        % Viscoelastic viscosity = (1-xelvis)*ETA
        etas0 = G.etas1.*(1-xelvis);
        % Vsicoelastic stress = xelvis*Sxy
        sxy0 = G.sxy1.*xelvis;
        
        % Normal stress
        %Viscoelasticity factor
        xelvis = G.etan1./(G.etan1 + timestep*G.mun1);
        % Viscoelastic viscosity = (1-xelvis)*ETA
        etan0 =G.etan1.*(1-xelvis);
        % Vsicoelastic stress = xelvis*Sxx
        sxx0 = G.sxx1.*xelvis;
    else
        etas0 = G.etas1;
        etan0 = G.etan1;
        % Stress will be calculated after Stokes solution
    end
    
    % Computing right part of mechanical viscoelastic equations
    % x-Stokes
    RX1=zeros(ynum+1,xnum);
    % y-Stokes
    RY1=zeros(ynum,xnum+1);
    % continuity
    RC1=zeros(ynum-1,xnum-1);
    
    % Grid points cycle
    for i=2:ynum
        for j=2:xnum
            
            % Right part of x-Stokes Equation
            if(j<xnum)
                RX1(i,j)=-gx*(G.rho1(i,j)+G.rho1(i-1,j))/2;  
                if viscoelastic
                    % Adding xelvis*dSxx0/dx
                    RX1(i,j)=RX1(i,j)-(sxx0(i-1,j)-sxx0(i-1,j-1))/xstpc1(j);
                    % Adding xelvis*dSxy0/dy
                    RX1(i,j)=RX1(i,j)-(sxy0(i,j)-sxy0(i-1,j))/ystp1(i-1);
                end
            end
            
            % Right part of y-Stokes Equation
            if(i<ynum)
                RY1(i,j)=-gy*(G.rho1(i,j)+G.rho1(i,j-1))/2;
                if viscoelastic
                    % Adding xelvis*dSyy0/dy using that Syy0=-Sxx0 (deviatoric stress)
                    RY1(i,j)=RY1(i,j)+(sxx0(i,j-1)-sxx0(i-1,j-1))/ystpc1(i);
                    % Adding xelvis*dSyx0/dx using that Syx0=Sxy0
                    RY1(i,j)=RY1(i,j)-(sxy0(i,j)-sxy0(i,j-1))/xstp1(j-1);
                end
            end
        end
    end
    
    comp_time(5) = toc;
    
    % Computing velocity field
    if (movemod==0)
        % Solving of Stokes and Continuity equations on nodes
        % and computing residuals
        % by calling function Stokes_Continuity_solver_grid()
        % with viscoelastic numerical viscosity
        % and modified right parts       
        [G.vx1,resx1,G.vy1,resy1,G.pr1,resc1]=newFast_Stokes_Continuity_solver_sandbox2(prfirst,etas0,etan0,xnum,ynum,gridx,gridy,RX1,RY1,RC1,bleft,bright,btop,bbottom,bintern,useSuiteSparse);
       % [G.vx1,resx1,G.vy1,resy1,G.pr1,resc1]=Fast_Stokes_Continuity_solver_sandbox(prfirst,etas0,etan0,xnum,ynum,gridx,gridy,RX1,RY1,RC1,bleft,bright,btop,bbottom,bintern);
        disp(['Stokes Relative Residuals: x, y, cont: ',num2str([max(resx1(:)) max(resy1(:)) max(resc1(:))])])
        
    end
    
    comp_time(6) = toc;
    
    % Solid body rotation
    if (movemod==1)
        for i=1:ynum+1
            for j=1:xnum+1
                % Vx
                if(j<xnum+1)
                    % Relative distance of vx node from the model center
                    dx=((j-1)*xstp1(j)-xsize/2)/(xsize/2);
                    dy=((i-1.5)*ystp1(i)-ysize/2)/(xsize/2);
                    dr=(dx^2+dy^2)^0.5;
                    % Set vx
                    G.vx1(i,j)=-vyright*dy;
                end
                % Vy
                if(i<ynum+1)
                    % Relative distance of vy node from the model center
                    dx=((j-1.5)*xstp1(j)-xsize/2)/(xsize/2);
                    dy=((i-1)*ystp1(i)-ysize/2)/(xsize/2);
                    dr=(dx^2+dy^2)^0.5;
                    % Set vy
                    G.vy1(i,j)=vyright*dx;
                end
            end
        end
    end    
    
    % Computing EPS'xx=-EPS'yy, EPSxy=EPSyx deviatoric strain rate tensor components from vx, vy
    
    % Computing spin tensor Espin
    G.exy = zeros(ynum,xnum);
    G.exx = zeros(ynum-1,xnum-1);
    G.esp = zeros(ynum,xnum);
    G.eii = zeros(ynum-1,xnum-1);
    
    % Grid points cycle

    % EPS'xx=-EPS'yy=1/2(dvx/dx-dvy/dy)
    xs1 = repmat(xstp1',ynum-1,1);
    ys1 = repmat(ystp1,1,xnum-1);
    G.exx = 0.5*((G.vx1(2:ynum,2:xnum)-G.vx1(2:ynum,1:xnum-1))./xs1 - (G.vy1(2:ynum,2:xnum)-G.vy1(1:ynum-1,2:xnum))./ys1);
    
    xs1 = repmat(xstpc1',ynum,1);
    ys1 = repmat(ystpc1,1,xnum);
    % EPSxy=EPSyx=1/2(dvx/dy+dvy/dx)
    % Espin=1/2(dvy/dx-dvx/dy) i.e. positive for clockwise rotation
    % (when x axis is directed rightward and y axis is directed downward)
    G.exy = 0.5*((G.vx1(2:ynum+1,1:xnum)-G.vx1(1:ynum,1:xnum))./ys1 + (G.vy1(1:ynum,2:xnum+1)-G.vy1(1:ynum,1:xnum))./xs1);
    G.esp = 0.5*((G.vy1(1:ynum,2:xnum+1)-G.vy1(1:ynum,1:xnum))./xs1 - (G.vx1(2:ynum+1,1:xnum)-G.vx1(1:ynum,1:xnum))./ys1);   
     
    % EPSii=(EPS'xx^2+EPSxy^2)^0.5
    G.eii = sqrt(G.exx(1:ynum-1,1:xnum-1).^2 + (G.exy(1:ynum-1,1:xnum-1).^2 + G.exy(2:ynum,1:xnum-1).^2 + G.exy(1:ynum-1,2:xnum).^2 + G.exy(2:ynum,2:xnum).^2)/4);
    
    % Check maximal velocity, ignoring the sticky air/water areas
    k = G.he>0.5;
    k = [G.he>0.5;k(1,:)];
    vxmax = max(abs(G.vx1.*k));
    k = G.he>0.5;
    k = [k(:,1) G.he>0.5];
    vymax = max(abs(G.vy1.*k),[],2);
    vxmax = interp1(gridx,vxmax',gridcx(2:end-1));
    vymax = interp1(gridy,vymax, gridcy(2:end-1));
    
    txmin = min(markmax*xstp1./vxmax);
    tymin = min(markmax*ystp1./vymax);
    
    disp(['Max velocity: ',num2str(max([vxmax(:);vymax(:)])*100*yr2sec),' cm/yr'])
    
    % The displacement time stepa
    timestep = min([txmin tymin Par.timemax timestep]);
    
    disp(['Displacement time step: ',num2str(timestep/yr2sec/1000),' ky'])
    
    if ntimestep == 137
        pause(1)
    end
    if viscoelastic
        
        % Computing new stresses and stress change using the displacement timestep
   
        % Shear stress
        % Viscoelasticity factor
        xelvis = G.etas1./(G.etas1 + timestep*G.mus1);
        % New viscoelastic stress = (1-xelvis)*2*ETA*EPSxy + xelvis*Sxy0
        sxy2  = 2*(1-xelvis).*G.etas1.*G.exy + xelvis.*G.sxy1;
        
        % Normal stress
        %Viscoelasticity factor
        xelvis = G.etan1./(G.etan1 + timestep*G.mun1);
        % New viscoelastic stress = (1-xelvis)*2*ETA*EPSxx + xelvis*Sxx0
        sxx2 = 2*(1-xelvis).*G.etan1.*G.exx + xelvis.*G.sxx1;

        % Elastic Stress change
        dsxy = sxy2-G.sxy1;
        dsxx = sxx2-G.sxx1;
        
        comp_time(7) = toc;
        
        % Computing strain rate and pressure for markers
        [MPR,MRAT,MEXX,MEXY] = mexmarker_pressure_strainrate(MCXN,MCYN,MCDX,MCDY,MXN,MYN,MDX,MDY,META,MMU(MI),sxx2,dsxx,sxy2,dsxy,G.exx,G.exy,G.pr1,timestep);
        
        % Computing subgrid stress changes for markers
        [MSXX,MSXY] = mexsubgridS_diffusion(MXN,MYN,MCXN,MCYN,MDX,MDY,MCDX,MCDY,META,MMU(MI),MSXX,MSXY,G.sxx1,G.sxy1,dsxx,dsxy,dsubgrids,timestep);
        
        comp_time(8) = toc;
        
    else
        
        comp_time(7) = toc;
        
        % Compute stress  S'xx
        sxx2 = 2*G.etan1.*G.exx;
        G.sxx1 = sxx2;
        
        % Compute shear stress S'xy
        sxy2 = 2*G.etas1.*G.exy;
        G.sxy1 = sxy2;
        
        % Computing strain rate and pressure for markers
        [MPR,MSXX,MSXY,MEXX,MEXY] = mexmarker_pressure_stress_strainrate(MCXN,MCYN,MCDX,MCDY,MXN,MYN,MDX,MDY,G.sxx1,G.sxy1,G.exx,G.exy,G.pr1);
        
       comp_time(8) = toc;
        
    end
    
    
    % Solving Temperature equation
    temp_problem = false;

    
    if (timestep>0 && tempmax>0)
        
        % Computing right part of temperature equation, radiogenic heat and
        % latent heat from phase change
        
        RT1 = G.hr1;        
        
        % Adiabatic heating on(1)/off(0)
        % Adding alp*T*DP/dt where DP/dt ~ vx*gx*rho+vy*gy*rho
        if (adiabyn==1)
            
            RT1(2:end-1,2:end-1) = RT1(2:end-1,2:end-1) + 0.5*G.ha1(2:end-1,2:end-1).*G.rho1(2:end-1,2:end-1).*(gx*(G.vx1(2:end-2,2:end-1) + G.vx1(3:end-1,2:end-1)) + gy*(G.vy1(2:end-1,2:end-2) + G.vy1(2:end-1,3:end-1)));
           
        end
        
        % Computing viscoelastic shear heating for Temperature nodes
        % Hs=2*Sxx*Sxx/2/etan+2*Sxy*Sxy/2/etas
        % Shear heating on(1)/off(0)
        % Limit to interior of the grid, no shear
        % heating 5 grid points in from left and right boundaries
        % and 5 grid points from bottom
        if(frictyn==1)
                        
            sh1 = sxx2.^2 ./G.etan1;
            RT1(6:ynum-5,6:xnum-5) = RT1(6:ynum-5,6:xnum-5) + sxy2(6:ynum-5,6:xnum-5).^2./G.etas1(6:ynum-5,6:xnum-5) + (sh1(5:ynum-6,5:xnum-6) + sh1(6:ynum-5,5:xnum-6) + sh1(5:ynum-6,6:xnum-5) + sh1(6:ynum-5,6:xnum-5))/4;
                        
        end
        
        % Solving temperature equation making (if needed) several thermal
        % timesteps for one displacement timestep
        % Set current thermal timestep
        timestept=timestep;
        % Set total thermal timestep
        timesteps=0;
        % Set old Temperature
        tk0=G.tk1;
        
        % Control break
        min_timestept = 365*24*3600; % 1 yr
        
        while (timesteps<timestep)
            
            % Solving Temperature equation with thermal timestep
            
            [tk2,rest]=Fast_Temperature_solver_grid(timestept,xnum,ynum,gridx,gridy,G.kt1,G.rhocp1,tk0,RT1,bleftt,brightt,btopt,bbottomt,useSuiteSparse);
            
            % Computing temperature changes
            dtk1=tk2-tk0;
            % Checking temperature changes
            dtkmax=max(max(abs(dtk1)));
            % Repeating temperature solution if temperature changes are too big
            
            if (dtkmax>abstempmax)
                disp(['Max. Absolute Temperature Change Exceeded: ',num2str(dtkmax),' > ',num2str(abstempmax), 'K'])
                disp('Unphysical temperature reached, no chance of convergence')
                temp_problem = true;
                break
            end
            
            if(dtkmax>tempmax)
                % Computing reduced timestep
                timestept = timestept*tempmax/dtkmax;
                disp(['Max. Temperature Change Exceeded: ',num2str(dtkmax),' > ',num2str(tempmax), 'K'])
                disp(['Temperature reduced timestep: ',num2str(timestept/yr2sec/1000),' ky'])
                % Solving Temperature equation with reduced timestep
                [tk2,rest]=Fast_Temperature_solver_grid(timestept,xnum,ynum,gridx,gridy,G.kt1,G.rhocp1,tk0,RT1,bleftt,brightt,btopt,bbottomt,useSuiteSparse);
                % Computing temperature changes
            end
            % Add total thermal timestep
            timesteps = timesteps + timestept;
            % Compute current thermal timestep
            if (timestept>timestep-timesteps)
                timestept=timestep-timesteps;
            end
            % Update old temperature
            tk0=tk2;
        end
        
        if temp_problem
            disp('Temperature solution failed to converge')
            break
        end
        
        % Compute temperature changes
        dtk1=tk2-G.tk1;
        
        % Computing subgrid diffusion for markers
       MTK0 = MTK;              
       MTK = mexsubgridT_diffusion(MXN,MYN,MDX,MDY,MTK0,G.tk1,G.rhocp1,dtk1,G.kt1,xstp1,ystp1,dsubgridt,timestep);
            
    end
    
    comp_time(9) = toc;
    
    % Moving Markers by velocity field
    if(markmove>0)
        % Create arrays for velocity and spin of markers
                 
         [MX,MY,MXN,MYN,MCXN,MCYN,MSXX,MSXY] = ...
             mexmovemarkersFast(MX,MY,MXN,MYN,MCXN,MCYN,MSXX,MSXY,...
             G.vx1,G.vy1,G.esp,gridx,gridy,gridcx,gridcy,markmove,timestep,xstp1,ystp1,xstpc1,ystpc1);
                  
         % Adding marker plastic strain based on grid strain rates
         meii = timestep*sqrt(MEXX.*MEXX + MEXY.*MEXY);
         k = plastic_yield > 0;
         MGII(k) = MGII(k) + meii(k);
        
         % Apply plastic strain healing by reducing the accumulated plastic
         % strain at every time step
         MGII = MGII/(1+timestep/tau_plastic);

         % Adding marker bulk strain based on grid strain rates        
         MBII = MBII + meii;
    end
    
    % Solide state phase transformations
    if  RxType
        
        % Eclogite phase
        if ntimestep > 3
            MECL = eclogite(RxType,MTK,MPR,MI,MECL,MPHASE,RGAS,timestep);
        end
        
        % 410 and 660 phases and rock-type change
        [MPH410,MPH660,MI] = MantlePhaseChange(MPHASE,MTK,MPR,MI,MPH410,MPH660,MRHOCUR,timestep,phase_function,MY);      
        
        % Zero out the latent heat effect of the initial distribution
        % We don't want a burst of heat from all the phase transistions,
        % i.e. we assume that transizition zone and lower mantle have
        % already transformed, it takes about 3 time steps to stabalize the
        % depth of the 660 and about 2 time steps for the 410.
        if ntimestep < 4
            MPH410(:,2) = 0;
            MPH660(:,2) = 0;
            
        end
       
    end
        
    if sticky_layer > 0       

        % Recomputing topography surface
        if strcmpi(LEMpar.apply_surfacepro,'fastscape')
            [gridt,H,E,LEMpar] = surface_topography_FS(gridx,gridy,xstp1,ystp1,gridcx,gridcy,xstpc1,ystpc1,G.vx1,G.vy1,gridt,H,E,tstp,LEMpar,timestep);
        else
            [gridt,H,E] = surface_topography3(gridx,gridy,xstp1,ystp1,gridcx,gridcy,xstpc1,ystpc1,G.vx1,G.vy1,gridt,H,E,tstp,LEMpar,timestep);
        end

        if any(isnan(gridt(:)))
           warning('NaN found in topography array')
        end
        
        % Take care of "crossing over"
        for mm1 = 1:marknum
            
            % Erosion-sedimentation
            % Find topography node for the marker
            xn=double(int16((MX(mm1)-gridt(1,1))/tstp-0.5))+1;
            if (xn<1)
                xn=1;
            end
            if (xn>tnum-1)
                xn=tnum-1;
            end
            % Save horizontal index
            TXN(mm1)=xn;
            
        end
        
        % Convert to sediment markers, using basement elevation
        
        % Compute relative distance to topography node
        dx = (MX-gridt(1,TXN)')/tstp;
        % Compute basement elevation above the marker
        dyB = H(1,TXN)'.*(1-dx)+H(1,TXN+1)'.*dx;
 
        % Compute topograhy elevation above the marker
        dyT = gridt(2,TXN)'.*(1-dx)+gridt(2,TXN+1)'.*dx;

        % water/air/rock that is not already sediment transform to sediment
        k = MY < dyB & MI ~=3 & MY > dyT; %(MI == 1 | MI == 2) & MY > dy;
        MI(k) = 3;   % Change marker type to sediment
        MRAT(k) = 1; % Reset strain rate ratio
        MGII(k) = 0; % Reset plastic strain
        MBII(k) = 0; % Reset bulk strain

%         % Compute topograhy elevation above the marker
%         dy = gridt(2,TXN)'.*(1-dx)+gridt(2,TXN+1)'.*dx;
% 
%         % Rocks to air transformation
%         k = MI > 1 & MI~=2 & MY < dy;
%        
%         MI(k) = 1;   % Change marker type
%         MRAT(k) = 1; % Reset strain rate ratio
%         MGII(k) = 0; % Reset plastic strain
%         MBII(k) = 0; % Reser bulk strain
        

        % Used for simulations of contintental basins containing lakes,
        % ensuring a water_level
        if ~isempty(LEMpar.sea_level)
            water_lev = max(gridt(2,:)) - LEMpar.max_depth;
            water_lev = max(LEMpar.offset-LEMPar.max_sea_level,water_lev);
            LEMpar.sea_level = LEMpar.offset-water_lev;
        end
        
        % Depricate function, kept for backward compatibility
        if water_depth>0 && timesum/yr2sec/1e6 > 1
            
            water_lev = get_water_level(gridt,water_depth,sticky_layer);
            
        end
        
        % Air to water transformation and vice versa (can't uplift water)

        if ~isempty(water_lev)
            k = MI == 1 & MY > water_lev;
            MI(k) = 2;
            MRAT(k) = 1; % Reset strain rate ratio
            MGII(k) = 0; % Reset plastic strain
            MBII(k) = 0; % Reset bulk strain
            k = MI == 2 & MY < water_lev;
            MI(k) = 1;
            MRAT(k) = 1; % Reset strain rate ratio
            MGII(k) = 0; % Reset plastic strain
            MBII(k) = 0; % Reset bulk strain
        end
        
        % Reset temperature to ttop for all sticky_markers, 
        
        % !!!!!!!!!!!!
        % should improve temperature stability to comment out -- needs to be further
        % investigated
        % !!!!!!!!!!!!        
        
        % k = MI == 2 | MI == 1;
        % MTK(k) = ttop;
        
        % Set the next horizon
        if timesum+timestep >= t_horizon
            i_horizon = i_horizon + 1;
            H(i_horizon,:) = gridt(2,:);
            t_horizon = t_horizon + dt_horizon;
        end
        
    end
    
    comp_time(10) = toc;

    % Save before adding markers
    % Also save MRHO and MLAMBDA
    if ntimestep==1  % Safe intial model
        
        time_outb = Par.yfreq_big;    % full size (big)
        time_outs = Par.yfreq_small;  % zoomed-in version (small)
        out_countb = 1;
        out_counts = 1;
        
        if RxType
        save([outpath,'/markers_',num2str(out_countb),'.mat'],...
            'MTK','MI','MX','MY','MRHO','MRHOCUR','MFLOW','MMU','MPL','MCP','MKT','MHR','MXN',...
            'MYN','MSXX','MSXY','META','MEXX','MEXY','MEXT','MEXTC','MPR','MGII','MBII','MRAT','MHE',...
            'MLAMBDA','MECL','MPH410','MPH660','marknum','timesum');
        else
            save([outpath,'/markers_',num2str(out_countb),'.mat'],...
            'MTK','MI','MX','MY','MRHO','MRHOCUR','MFLOW','MMU','MPL','MCP','MKT','MHR','MXN',...
            'MYN','MSXX','MSXY','META','MEXX','MEXY','MEXT','MEXTC','MPR','MGII','MBII','MRAT','MHE',...
            'MLAMBDA','marknum','timesum');
        end

        if strcmpi(LEMpar.apply_surfacepro,'fastscape')
            FastScape.topo = reshape(double(LEMpar.hfinal),LEMpar.nx,LEMpar.ny)';
            FastScape.basement = reshape(double(LEMpar.bfinal),LEMpar.nx,LEMpar.ny)';
            FastScape.lake_depth = reshape(double(LEMpar.lake_depth),LEMpar.nx,LEMpar.ny)';
            FastScape.nx = LEMpar.nx;
            FastScape.ny = LEMpar.ny;
            FastScape.xl = LEMpar.xl;
            FastScape.yl = LEMpar.yl;
        else
            FastScape = [];
        end

        save([outpath,'/grids_',num2str(out_countb),'.mat'],...
            'G','gridt','H','E','FastScape','timesum','water_lev','comp_time')
        
        if addmarkers
            save([outpath,'/bounds.mat'],'btop','bbottom','bleft','bright','btopt','bbottomt','bleftt','brightt','bintern','reservoir')
        else
            save([outpath,'/bounds.mat'],'btop','bbottom','bleft','bright','btopt','bbottomt','bleftt','brightt','bintern')
        end
        
        % small box
        k = MX >= Par.markx(1) & MX <= Par.markx(2) & MY >= Par.marky(1) & MY <= Par.marky(2);
        
        % Should MRHOCUR and MLAMBDA be saved here?
        save_markersbox(outpath,k,MTK,MI,MX,MY,MRHO,MFLOW,MMU,MPL,MCP,MKT,MHR,...
            MSXX,MSXY,META,MEXX,MEXY,MPR,MGII,MBII,MRAT,out_counts,timesum);
        
        % save  limits of small output box and topo grid of all the model:
        markx = Par.markx; marky= Par.marky;
        save([outpath,'/grids_',num2str(out_counts),'s.mat'],...
            'markx','marky','gridt','timesum');
    end
    
    
    % Advance in time
    timesum = timesum+timestep;
    
    % Save acording to step frequency (true) or time frequency (false)
    write_big = false;
    write_small = false;
    
    if Par.outtype  % Step
        
        if ~mod(ntimestep,Par.nfreq_big)
            write_big = true; 
            out_countb = ntimestep;        
        end
        if ~mod(ntimestep,Par.nfreq_small)
            write_small = true; 
            out_counts = ntimestep;         
        end
        
    else    % Time
        
        if  timesum >= time_outb
            time_outb = time_outb + Par.yfreq_big;
            write_big = true;
            out_countb = out_countb + 1;
        end
        
        if  timesum+timestep >= time_outs
            time_outs = time_outs + Par.yfreq_small;
            write_small = true;
            out_counts = out_counts + 1;
        end
       
    end
            
    if write_big    
        % Full Size

        if RxType
            save([outpath,'/markers_',num2str(out_countb),'.mat'],...
                'MTK','MI','MX','MY','MRHO','MRHOCUR','MFLOW','MMU','MPL','MCP','MKT','MHR','MXN',...
                'MYN','MSXX','MSXY','META','MEXX','MEXY','MEXT','MEXTC','MXM','MPR','MGII','MBII','MRAT',...
                'MHE','MLAMBDA','MECL','MPH410','MPH660','marknum','timesum');
        else
            save([outpath,'/markers_',num2str(out_countb),'.mat'],...
                'MTK','MI','MX','MY','MRHO','MRHOCUR','MFLOW','MMU','MPL','MCP','MKT','MHR','MXN',...
                'MYN','MSXX','MSXY','META','MEXX','MEXY','MEXT','MEXTC','MXM','MPR','MGII','MBII','MRAT',...
                'MHE','MLAMBDA','marknum','timesum');
        end

        % You can't save python modules
        if strcmpi(LEMpar.apply_surfacepro,'fastscape')
            FastScape.topo = reshape(double(LEMpar.hfinal),LEMpar.nx,LEMpar.ny)';
            FastScape.basement = reshape(double(LEMpar.bfinal),LEMpar.nx,LEMpar.ny)';
            FastScape.lake_depth = reshape(double(LEMpar.lake_depth),LEMpar.nx,LEMpar.ny)';            
            FastScape.nx = LEMpar.nx;
            FastScape.ny = LEMpar.ny;
            FastScape.xl = LEMpar.xl;
            FastScape.yl = LEMpar.yl;
        else
            FastScape = [];
        end

        save([outpath,'/grids_',num2str(out_countb),'.mat'],...
            'G','gridt','H','E','FastScape','timesum','water_lev','comp_time');
    end
    
    if  write_small
        % Zoomed-in 
        
        k = MX >= Par.markx(1) & MX <= Par.markx(2) & MY >= Par.marky(1) & MY <= Par.marky(2);
        
        % Should MRHOCUR and MLAMBDA be saved here?
        save_markersbox(outpath,k,MTK,MI,MX,MY,MRHO,MFLOW,MMU,MPL,MCP,MKT,MHR,...
            MSXX,MSXY,META,MEXX,MEXY,MPR,MGII,MBII,MRAT,out_counts,timesum);
        
        % save  limits of small output box and topo grid of all the model:
        save([outpath,'/grids_',num2str(out_counts),'s.mat'],...
            'markx','marky','gridt','timesum');
        
    end
    
    % If move inner boundary condition
    if movebound
        % position of interior boundary
        ix = ix + con_velo*timestep;
        disp(['Innerbound Position: ',num2str(ix/1000)])
        
        % Find if we have moved passed the next "node"
        k = bintern(1);
        if ix >= G.gridx(k+1)
            k = k + 1;
            bintern(1)=k;
            disp('Innerbound Position Updated')
        end
        
    end
    
    comp_time(11) = toc;

    % Fill Empty Cell
    mx = []; my = []; mi = []; mxn = []; myn = []; newmarkers = 0;
    
    if markerreflect && ~fillemptycells && ~addmarkers
        
        % Take outflow markers and reflect across as inflow markers
        [mx,my,mxn,myn,mcxn,mcyn,mi,meta,mtk] = reflectmarkers(gridx,gridy,gridcx,gridcy,MX,MY,MI,META,MTK);        
        newmarkers = length(mx);
        
    end
    
    if fillemptycells && ~markerreflect && ~addmarkers
                
        [mx,my,mi,mxn,myn,mcxn,mcyn,newmarkers] = fillemptycell(gridx,gridy,gridcx,gridcy,MXN,MYN,mxcell,mycell,filloption,newmarkertype);
        
    end
    
    % Inject markers at left and right boundaries
    if addmarkers && ~markerreflect && fillemptycells
        
        bxl = bxl + vxleft*timestep;
        bxr = bxr + vxright*timestep;
        
        if bxl > xstp1(1)/8
            bl = true;
            bxl = 0;
        else
            bl = false;
        end
        
        if bxr > xstp1(end)/8
            br = true;
            bxr = 0;
        else
            br = false;
        end
        
        [mx,my,mxn,myn,mcxn,mcyn,mi] = injectmarkers(gridx,gridy,gridcy,BCL,BCR,mx,my,mxn,myn,mcxn,mcyn,mi,reservoir,bl,br);
        newmarkers = length(mi);
        mx = mx(:); my = my(:); mxn = mxn(:);myn = myn(:);
        mcxn = mcxn(:); mcyn = mcyn(:);
        mi=mi(:);
        
    end
    
    % Compute pressure, rheology, stress, strain and temperature of new
    % markers
    
    if newmarkers > 0

        zmarks = zeros(size(mx));
        msxx = zmarks;
        msxy = zmarks;
        
        % Define normalized distances from new markers to the upper left node;
        mdx = (mx - gridx(mxn))./xstp1(mxn);
        mdy = (my - gridy(myn))./ystp1(myn);
        mcdx = (mx - gridcx(mcxn+1))./xstpc1(mcxn+1);
        mcdy = (my - gridcy(mcyn+1))./ystpc1(mcyn+1);
        
        % Interpolate temperature and viscosity to new markers
        if ~markerreflect
            mtk = zmarks;
            meta = zmarks;

            for mm1 = 1:newmarkers

                xn = mxn(mm1);
                yn = myn(mm1);

                % Define normalized distances from marker to the upper left node;
                dx = mdx(mm1);
                dy = mdy(mm1);

                % Calculate Marker temperature from four surrounding nodes
                tkm = (1.0-dx).*(1.0-dy).*G.tk1(yn,xn);
                tkm = tkm + (1.0-dx).*dy.*G.tk1(yn+1,xn);
                tkm = tkm + dx.*(1.0-dy).*G.tk1(yn,xn+1);
                tkm = tkm + dx.*dy.*G.tk1(yn+1,xn+1);
                mtk(mm1) = tkm;

                % Calculate Marker viscosity from four surrounding nodes
                % Using arithmetic mean -- could use geomteric instead?
                eta = (1.0-dx).*(1.0-dy).*G.etas1(yn,xn);
                eta = eta + (1.0-dx).*dy.*G.etas1(yn+1,xn);
                eta = eta + dx.*(1.0-dy).*G.etas1(yn,xn+1);
                eta = eta + dx.*dy.*G.etas1(yn+1,xn+1);
                meta(mm1) = eta;
            end

            % Comput subgrid temperature diffusion
            mtk0 = mtk;
            mtk = subgridT_diffusion(mxn,myn,mdx,mdy,mtk0,G.tk1,G.rhocp1,dtk1,G.kt1,xstp1,ystp1,dsubgridt,timestep);

        end

        if viscoelastic

            % Computing strain rate and pressure for new markers
            [mpr,mrat,mexx,mexy] = mexmarker_pressure_strainrate(mcxn,mcyn,mcdx,mcdy,mxn,myn,mdx,mdy,meta,MMU(mi),sxx2,dsxx,sxy2,dsxy,G.exx,G.exy,G.pr1,timestep);
            
            % Computing subgrid stress changes for markers
            [msxx,msxy] = mexsubgridS_diffusion(mxn,myn,mcxn,mcyn,mdx,mdy,mcdx,mcdy,meta,MMU(mi),msxx,msxy,G.sxx1,G.sxy1,dsxx,dsxy,dsubgrids,timestep);
            
        else
            
            % Computing strain rate and pressure for markers
             mrat = zmarks + 1;
            [mpr,msxx,msxy,mexx,mexy] = mexmarker_pressure_stress_strainrate(mcxn,mcyn,mcdx,mcdy,mxn,myn,mdx,mdy,G.sxx1,G.sxy1,G.exx,G.exy,G.pr1);
            
        end
                       
        MX = [MX; mx];
        MY = [MY; my];
        MSXX = [MSXX; msxx];
        MSXY = [MSXY; msxy];
        MTK = [MTK; mtk];
        MI = [MI; mi];
        
        MXN = [MXN; mxn];
        MYN = [MYN; myn];
        MCXN = [MCXN; mcxn];
        MCYN = [MCYN; mcyn];
        MRAT = [MRAT; 1+zmarks];
        MGII = [MGII; zmarks];
        MBII = [MBII; zmarks];
        MPR = [MPR; mpr];
        META = [META; meta];
        MEXX = [MEXX; mexx];
        MEXY = [MEXY; mexy];
        MXM = [MXM; zmarks];
        MEXT = [MEXT; zmarks];
        MEXTC = [MEXTC; zmarks];
        TXN = [TXN zmarks'];
        MLAMBDA = [MLAMBDA; 1+zmarks];
        if RxType
            MECL = [MECL; zmarks zmarks+1e4 zmarks zmarks zmarks];
            MPH410 = [MPH410; zmarks zmarks];
            MPH660 = [MPH660; zmarks zmarks];
        end
        plastic_yield = [plastic_yield; zmarks];
        marknum = marknum + newmarkers;
        markers_added = markers_added + newmarkers;
    end    
    
    
    disp(['Model time: ',num2str(timesum/yr2sec/1000),' ky'])
    
    if timesum >= Par.modeltime
        break
    end 
    
    comp_time(12) = toc;
    
    % Reset the restart
    restart = false;
    
    % If the grid is being stretched with no regridding, increase the
    % topography tracking surface as needed and adjust the boundary
    % conditions (top and bottom) to conserve mass and to make sure that no
    % adiabatic heating occurs as a result of moving the bottom thermal
    % boundary upwards (or downwards).

    if moving_grid && ~regridx
 
        % Defining new gridline positions for irregular basic grid
        % Computing model size change
        xsize0 = G.xsize;
        ysize0 = G.ysize;
        gridx1 = G.gridx(1);  % Current x(1) coordinate
        
        % Increase/decrease x-dimension of grid
        dxsize = (vxright-vxleft)*timestep;
        xsize = xsize0 + dxsize; % New Horizontal size
        
        % Decrease/increase y-dimension of grid (to keep area constant)
        ysize = (xsize0/xsize)*ysize0; % Vertical size
        
        
        G.xsize = xsize;
        G.ysize = ysize;
        tmpG = construct_grid(G);
        gridx = tmpG.x'; gridy = tmpG.y';
        clear tmpG

        % Ensure we have the correct grid size in real-space
        G.xsize = gridx(end);
        G.ysize = gridy(end);
        G.gridy = gridy;

        % Now recenter the grid by shifting to the left (note that if
        % vxleft .NEQ. vxright, the high-resolution won't cover the same
        % area, and will drift over time favoring the lower velocity side)
        G.gridx = gridx  + vxleft*timestep + gridx1;
        gridx1 = G.gridx(1);
        gridx1N = G.gridx(end);
        gridx = G.gridx;
      
        % Add a topography node when the new grid extends a distance past the
        % topography grid resolution; topography grid should always
        % encompas the model grid (assumes the grid will not move more than
        % tstp per step).

        % Not sure if we also need to remove topography nodes if
        % grid is compressed
        adjust_gridt = [false false];
        if (gridt(1,1)-gridx1)>0   % Left of model space
            gridt = [gridt(:,1) gridt];
            gridt(1,1) = gridt(1,1) - tstp;
            tnum = tnum + 1;
            H = [H(:,1) H];
            E = [E(1) E];
            adjust_gridt(1) = true;
        end
        if (gridx1-gridt(1,end)) > 0   % Right of model space
            gridt = [gridt gridt(:,end)];
            gridt(1,end) = gridt(1,end) + tstp;
            tnum = tnum + 1;
            H = [H H(:,end)];
            E = [E E(end)];
            adjust_gridt(2) = true;
        end
 
        % Reinitialize FastScape
        if strcmpi(LEMpar.apply_surfacepro,'fastscape') && any(adjust_gridt)
            hinit = reshape(double(LEMpar.hfinal),LEMpar.nx,LEMpar.ny)';
            binit = reshape(double(LEMpar.bfinal),LEMpar.nx,LEMpar.ny)';

            if adjust_gridt(1)
                LEMpar.nx = LEMpar.nx + 1;
                hinit = [hinit(:,1) hinit];
                binit = [binit(:,1) binit];
                
            end
            if adjust_gridt(2)
                LEMpar.nx = LEMpar.nx + 1;
                hinit = [hinit hinit(:,end)];
                binit = [binit binit(:,end)];
            end

            LEMpar.xl = abs(gridt(1,end) - gridt(1,1));
            LEMpar.nn = LEMpar.nx*LEMpar.ny;

            hinit = reshape(hinit',1,LEMpar.nn);
            binit = reshape(binit',1,LEMpar.nn);

            % Current step
            istep = double(LEMpar.FS.fastscape_get_step());

            % Re-initilize fastscape
            LEMpar.FS.fastscape_destroy()
            LEMpar.FS.fastscape_init()

            LEMpar.FS.fastscape_set_nx_ny(LEMpar.nx,LEMpar.ny);
            LEMpar.FS.fastscape_setup();
            LEMpar.FS.fastscape_set_xl_yl(LEMpar.xl,LEMpar.yl);

            % Surfrace processes parameters
            Kf = LEMpar.Kf + zeros(1,LEMpar.nn);
            Kd = LEMpar.Kd + zeros(1,LEMpar.nn);
            Kf = LEMpar.NP.array(Kf);
            Kd = LEMpar.NP.array(Kd);
            LEMpar.FS.fastscape_set_erosional_parameters(Kf,LEMpar.kfsed,LEMpar.m,LEMpar.n,Kd,LEMpar.kdsed,LEMpar.Fd,LEMpar.Fd,LEMpar.expp);

            LEMpar.FS.fastscape_set_marine_parameters( LEMpar.sea_level, LEMpar.p1, LEMpar.p2, LEMpar.z1, LEMpar.z2, LEMpar.r, LEMpar.L, LEMpar.kds1, LEMpar.kds2)

            % Set the topography
            LEMpar.hfinal = LEMpar.NP.array(hinit);
            LEMpar.FS.fastscape_init_h(LEMpar.hfinal);

            % Set the basement
            LEMpar.bfinal = LEMpar.NP.array(binit);
            LEMpar.FS.fastscape_set_basement(LEMpar.bfinal);

            % Initialize lake depth container
            LEMpar.lake_depth = LEMpar.NP.array(zeros(1,LEMpar.nn));

            % Precipitation
            p = LEMpar.p*ones(1,LEMpar.nn);
            p = LEMpar.NP.array(p);
            LEMpar.FS.fastscape_set_precip(p);

            % Set boundary conditions
            bc = 101; % Left and Right fixed, Top and Bottom cyclic
            LEMpar.FS.fastscape_set_bc(int32(bc));

        end
        % New Upper boundary: Free slip
        % vx(1,j)=btop(j,1)+vx(2,j)*btop(j,2)
        btop(:,1)=0;
        btop(:,2)=1;
        % vy(1,j)=btop(j,3)+vy(2,j)*btop(j,4)
        btop(:,3)= sticky_layer*(vxright-vxleft)/xsize;
        btop(:,4)= 0;
        
        % New Lower boundary: Free slip
        % vx(ynum+1,j)=bbottom(j,1)+vx(ynum,j)*bbottom(j,2)
        bbottom(:,1)= 0;
        bbottom(:,2)= 1;
        % vy(ynum,j)=bbottom(j,3)+vy(ynum-1,j)*bbottom(j,4)
        bbottom(:,3)= (vxright-vxleft)*(sticky_layer-ysize)/xsize;
        bbottom(:,4)= 0;
        
        % New thermal lower boundary conditions
        % Adiabatic T increase (0.5 K/km) in responce to model thickening
        % Temperature at bottom of model space
        
        if update_tbottom
            tbottom = tpotential + tgrad*(ysize-sticky_layer);
        end
        
        % Lower boundary
        % tk(ynum,j)=bbottomt(j,1)+tk(ynum-1,j)*bbottomt(j,2)
        bbottomt(:,1)= tbottom;
        bbottomt(:,2)= 0;
        
    end
    
    % Move the high resolution area
    if regridx && ~moving_grid
       
        
    end
        
    % Add the imposed horizontal velocity after some time in (in Myr)
    if add_time
        T_Ma = timesum/yr2sec/1e6;
        
        if T_Ma > add_time
            vxright = v1;
            vxleft = -v1;
            vxleft = vxleft/100/yr2sec;
            vxright = vxright/100/yr2sec;

            btop(:,3) = sticky_layer*(vxright-vxleft)/xsize;
            bbottom(:,3) = (vxright-vxleft)*(sticky_layer-ysize)/xsize;
            bleft(:,1) = vxleft;
            bright(:,1) = vxright;
            
            moving_grid = true;
            update_tbottom = true;            
            gridx0 = G.gridx(1);
            gridxN = G.gridx(end);
            
            add_time = false;  % no longer needed
        end
    end  
    % Remove the imposed horizontal velocity after some time (in Myr)
    if remove_time
        T_Ma = timesum/yr2sec/1e6;
        
        if T_Ma > remove_time
            vxright = 0;
            vxleft = 0;
            
            btop(:,3) = 0;
            bbottom(:,3) = 0;
            bleft(:,1) = 0;
            bright(:,1) = 0;
            
            moving_grid = false;
            remove_time = false;  % Since we did this, we can turn this off       
        end
    end    
    
    %--------------------------------------------------------
    % Add a plume
    %--------------------------------------------------------
    T_Ma = timesum/yr2sec/1e6;

    
    if add_plume
        
        if T_Ma >= add_plume
            
            % Check that parameters exist
            if ~exist('xplume','var') || ~exist('yplume','var') || ~exist('rplumex','var') || ~exist('rplumey','var') || ~exist('Texcess','var')
                error('If using add_plume, you must define plume options in input model script')
            end
            
            dx = MX - xplume;                      
            dy = MY - min(yplume,ysize);
            
            sigx = rplumex/sqrt(2*log(2));
            sigy = rplumey/sqrt(2*log(2));
                  
            % Gaussian
            switch plume_type
                
                case   'gaussian'
                    
                    dT = Texcess*exp(-0.5*((dx/sigx).^2 + (dy/sigy).^2));
                    k = dT >= 10;
                    MTK = MTK + dT;
                    % Plume Markers to track
                    MI(k) = 14;
                    
                case 'spherical'
                    
                    rplume2 = rplumex^2;
                    r2 = dx.^2 + dy.^2;
                    % Plume Markers to track
                    k = r2 <= rplume2;
                    MI(k) = 14;
                    MTK(k) = MTK(k) + Texcess;
                    
                otherwise
                    
                    error('Unknown plume type')

            end
                   
            add_plume = false;  % Since we did this, we can turn this off
            
             Par.timemax = 25e3*yr2sec;
             Par.yfreq_big = 50e3*yr2sec;  % Change ouput time
             time_outb = timesum + Par.yfreq_big;
%             
        end
        
    end
    
    
    comp_time(13) = toc;
    
    disp(['Run time: ',num2str(comp_time(end)),' sec'])
    
    comp_time(:) = NaN;
    
end

disp(['Markers added: ',num2str(markers_added)])
disp(['Markers removed: ',num2str(markers_removed)])

% Save the final output
% Include MRHO and MLAMBDA
if RxType
save([outpath,'/markers_',num2str(out_countb+1),'.mat'],...
    'MTK','MI','MX','MY','MRHO','MRHOCUR','MFLOW','MMU','MPL','MCP','MKT','MHR','MXN',...
    'MYN','MSXX','MSXY','META','MEXX','MEXY','MPR','MEXT','MEXTC','MXM','MGII','MBII','MRAT',...
    'MHE','MLAMBDA','MECL','MPH410','MPH660','marknum','timesum');
else
save([outpath,'/markers_',num2str(out_countb+1),'.mat'],...
    'MTK','MI','MX','MY','MRHO','MRHOCUR','MFLOW','MMU','MPL','MCP','MKT','MHR','MXN',...
    'MYN','MSXX','MSXY','META','MEXX','MEXY','MPR','MEXT','MEXTC','MXM','MGII','MBII','MRAT',...
    'MHE','MLAMBDA','marknum','timesum');    
end

% You can't save python modules
if strcmpi(LEMpar.apply_surfacepro,'fastscape')
    FastScape.topo = reshape(double(LEMpar.hfinal),LEMpar.nx,LEMpar.ny)';
    FastScape.basement = reshape(double(LEMpar.bfinal),LEMpar.nx,LEMpar.ny)';
    FastScape.lake_depth = reshape(double(LEMpar.lake_depth),LEMpar.nx,LEMpar.ny)';
    FastScape.nx = LEMpar.nx;
    FastScape.ny = LEMpar.ny;
    FastScape.xl = LEMpar.xl;
    FastScape.yl = LEMpar.yl;
    LEMpar.FS.fastscape_destroy()
    LEMpar = rmfield(LEMpar,{'NP';'FS';'hfinal';'bfinal';'lake_depth'});
else
    FastScape = [];
end

save([outpath,'/grids_',num2str(out_countb+1),'.mat'],...
    'G','gridt','H','E','FastScape','timesum','LEMpar','water_lev','comp_time','bleft','bright','btop','bbottom','bintern');

if temp_problem
    save([outpath,'debug_',num2str(ntimestep),'.mat'])
end
