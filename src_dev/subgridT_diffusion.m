function MTK = subgridT_diffusion(MXN,MYN,MDX,MDY,MTK,tk1,rhocp1,dtk1,kt1,xstp1,ystp1,dsubgridt,timestep)

% Clear subgrid temperature changes for nodes
[ynum,xnum]=size(tk1);
marknum = length(MTK);
dtkn=zeros(ynum,xnum);
% Clear wights for basic nodes
wtnodes=zeros(ynum,xnum);

% Marker cycle
for mm1 = 1:1:marknum
        
    % Interpolating temperature changes from basic nodes
    %
    % Define indexes for upper left node in the cell where the marker is
    xn=MXN(mm1);
    yn=MYN(mm1);
    
    % Define normalized distances from marker to the upper left node;
    dx=MDX(mm1);
    dy=MDY(mm1);
    
    % Compute marker weight koefficient from cell dimensions
    % Number of markers in a cell is in invert proportion to the cell volume
    mwt=1;%/xstp1(xn)/ystp1(yn);
    
    if mm1 == 3145
        mm1;
    end
    % Interpolate old nodal temperature for the marker
    tkm=0;
    tkm=tkm+(1.0-dx).*(1.0-dy).*tk1(yn,xn);
    tkm=tkm+(1.0-dx).*dy.*tk1(yn+1,xn);
    tkm=tkm+dx.*(1.0-dy).*tk1(yn,xn+1);
    tkm=tkm+dx.*dy.*tk1(yn+1,xn+1);
    % Calculate Nodal-Marker subgrid temperature difference
    dtkm=tkm-MTK(mm1);
    % Compute nodal k and RHO*Cp for the marker
    % k
    ktm=0;
    ktm=ktm+(1.0-dx).*(1.0-dy).*kt1(yn,xn);
    ktm=ktm+(1.0-dx).*dy.*kt1(yn+1,xn);
    ktm=ktm+dx.*(1.0-dy).*kt1(yn,xn+1);
    ktm=ktm+dx.*dy.*kt1(yn+1,xn+1);
    % RHO*Cp
    rhocpm=0;
    rhocpm=rhocpm+(1.0-dx).*(1.0-dy).*rhocp1(yn,xn);
    rhocpm=rhocpm+(1.0-dx).*dy.*rhocp1(yn+1,xn);
    rhocpm=rhocpm+dx.*(1.0-dy).*rhocp1(yn,xn+1);
    rhocpm=rhocpm+dx.*dy.*rhocp1(yn+1,xn+1);
    
    % Compute local thermal diffusion timescale for the marker
    tdm=rhocpm/ktm/(2/xstp1(xn)^2+2/ystp1(yn)^2);
    
    % Computing subgrid diffusion
    sdif=-dsubgridt*timestep/tdm;
    if(sdif<-30)
        sdif=-30;
    end
    dtkm=dtkm*(1-exp(sdif));
    
     %Computing new temperature for the marker
    if 273 - (MTK(mm1)+dtkm) > 0.001
        mm1;
    end
    % Correcting old temperature for the marker
    MTK(mm1)=MTK(mm1)+dtkm;
    
    
    % Interpolating subgrid temperature changes to 4 nodes
    dtkn(yn,xn)=dtkn(yn,xn)+(1.0-dx).*(1.0-dy).*dtkm*mwt;
    wtnodes(yn,xn)=wtnodes(yn,xn)+(1.0-dx).*(1.0-dy)*mwt;
    
    dtkn(yn+1,xn)=dtkn(yn+1,xn)+(1.0-dx).*dy.*dtkm*mwt;
    wtnodes(yn+1,xn)=wtnodes(yn+1,xn)+(1.0-dx).*dy*mwt;
    
    dtkn(yn,xn+1)=dtkn(yn,xn+1)+dx.*(1.0-dy).*dtkm*mwt;
    wtnodes(yn,xn+1)=wtnodes(yn,xn+1)+dx.*(1.0-dy)*mwt;
    
    dtkn(yn+1,xn+1)=dtkn(yn+1,xn+1)+dx.*dy.*dtkm*mwt;
    wtnodes(yn+1,xn+1)=wtnodes(yn+1,xn+1)+dx.*dy*mwt;
    
end

% Computing subgrid diffusion for nodes
k = wtnodes~=0;
dtkn(k) = dtkn(k)./wtnodes(k);

% Subtracting subgrid diffusion part from nodal temperature changes
dtk1=dtk1-dtkn;

tk1 = dtk1+tk1;
% Updating temperature for markers
for mm1 = 1:1:marknum
    
    % Interpolating temperature changes from basic nodes
    %
    % Define indexes for upper left node in the cell where the marker is
    xn=MXN(mm1);
    yn=MYN(mm1);
    
    % Define normalized distances from marker to the upper left node;
    dx=MDX(mm1);
    dy=MDY(mm1);
    
    % Calculate Marker temperature change from four surrounding nodes
    dtkm=0;
    dtkm=dtkm+(1.0-dx).*(1.0-dy).*dtk1(yn,xn);
    dtkm=dtkm+(1.0-dx).*dy.*dtk1(yn+1,xn);
    dtkm=dtkm+dx.*(1.0-dy).*dtk1(yn,xn+1);
    dtkm=dtkm+dx.*dy.*dtk1(yn+1,xn+1);    
     
    %Computing new temperature for the marker
    MTK(mm1)=MTK(mm1)+dtkm;        
    
    % Correct the minimum value
    min_TK = min(min(tk1([yn yn+1],[xn xn+1])));
    if MTK(mm1) < min(min(tk1([yn yn+1],[xn xn+1])))
        MTK(mm1) = min_TK;
    end
end
