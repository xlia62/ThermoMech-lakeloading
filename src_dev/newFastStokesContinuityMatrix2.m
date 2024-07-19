% Function Stokes_Continuity_solver_sandbox()
% This function formulates and solves  
% Stokes and Continuity equations defined on 2D staggered irregularly spaced grid
% with specified resolution (xnum, ynum) and grid lines positions (gridx, gridy)
% given distribution of right parts for all equations (RX,RY,RC) on the grid 
% and given variable shear (etas) and normal (etan) viscosity distributions 
% pressure is normalized relative to given value (prnorm) in the first cell
%
% Velocity Boundary condition specified by bleft,bright,btop,bbottom,bintern 
% are implemented from ghost nodes 
% directly into Stokes and continuity equations
%
% Function returns solution for velocity and pressure (vx,vy,pr)
% and distribution of residuals (resx,resy,resc)
function [I,J,L,R] = newFastStokesContinuityMatrix2(prfirst,etas,etan,xnum,ynum,gridx,gridy,RX,RY,RC,bleft,bright,btop,bbottom,bintern,pscale)

% Staggered Grid for Multigrid
% 
%     vx       vx       vx    
%
% vy  +---vy---+---vy---+   vy
%     |        |        |
%     vx   P   vx   P   vx    
%     |        |        |
% vy  +---vy---+---vy---+   vy
%     |        |        |
%     vx   P   vx   P   vx    
%     |        |        |
% vy  +---vy---+---vy---+   vy
%
%     vx       vx       vx    
% 
% Lines show basic grid
% Basic (density) nodes are shown with +
% Ghost nodes shown outside the basic grid
% are used for boundary conditions

% Boundary conditions
% Pressure boundary condition
% Pressure in first cell
bpres=0;
prnorm=prfirst(2);
% Channel flow top->bottom
if (prfirst(1)==1)
    bpres=1;
    prnorm=prfirst(2);
end

% Velocity boundary conditions
btopx(:,1)=btop(:,1);
btopx(:,2)=btop(:,2);
btopy(:,1)=btop(:,3);
btopy(:,2)=btop(:,4);
bbottomx(:,1)=bbottom(:,1);
bbottomx(:,2)=bbottom(:,2);
bbottomy(:,1)=bbottom(:,3);
bbottomy(:,2)=bbottom(:,4);
bleftx(:,1)=bleft(:,1);
bleftx(:,2)=bleft(:,2);
blefty(:,1)=bleft(:,3);
blefty(:,2)=bleft(:,4);
brightx(:,1)=bright(:,1);
brightx(:,2)=bright(:,2);
brighty(:,1)=bright(:,3);
brighty(:,2)=bright(:,4);

% Prescribed internal velocity condition


% Computing grid steps for basic nodes
xstp=zeros(xnum-1,1);
ystp=zeros(ynum-1,1);
for i=1:1:xnum-1
    xstp(i)=gridx(i+1)-gridx(i);
end
for i=1:1:ynum-1
    ystp(i)=gridy(i+1)-gridy(i);
end
% Computing grid steps for vx and vy nodes
xstpc=zeros(xnum,1);
ystpc=zeros(ynum,1);
% First and last steps (for ghost nodes)
xstpc(1)=xstp(1);
xstpc(xnum)=xstp(xnum-1);
ystpc(1)=ystp(1);
ystpc(ynum)=ystp(ynum-1);
for i=2:1:xnum-1
    xstpc(i)=(gridx(i+1)-gridx(i-1))/2;
end
for i=2:1:ynum-1
    ystpc(i)=(gridy(i+1)-gridy(i-1))/2;
end


% Average x and y steps
xstpavr=(gridx(xnum)-gridx(1))/(xnum-1);
ystpavr=(gridy(ynum)-gridy(1))/(ynum-1);

% Horizontal shift index
ynum3=(ynum-1)*3;

% Creating matrix
nzeros = xnum*ynum*26;
I = zeros(nzeros,1);
J = zeros(nzeros,1);
L = zeros(nzeros,1);
%nzeros = xnum*ynum*25; % to watch growth

%L=sparse((xnum-1)*(ynum-1)*3,(xnum-1)*(ynum-1)*3);
R=zeros((xnum-1)*(ynum-1)*3,1);

% Solving of Stokes and continuity equations on nodes

% Assemble the sparse "Mass" matrix to minimize the number of logical
% operations we break things up into several loops

%------------------------------------------------------------------
%% Top Left Corner
%------------------------------------------------------------------
c= 0;

j = 1;
i = 1;
ivx = ((j-1)*(ynum-1)+(i-1))*3+1;
ivy = ivx + 1;
ipr = ivx + 2;

% x-Stokes coefficients

% Right hand side
R(ivx,1) = RX(i+1,j+1);

% Central Vx node
c = c + 1;
I(c) = ivx;
J(c) = ivx;
L(c) = -2*(etan(i,j+1)/xstp(j+1)+etan(i,j)/xstp(j))/xstpc(j+1) ...
    -(etas(i+1,j+1)/ystpc(i+1)+etas(i,j+1)/ystpc(i))/ystp(i);

% Left Vx boundary condition
L(c) = L(c) + bleftx(i+1,2)*2*etan(i,j)/xstp(j)/xstpc(j+1);
R(ivx,1) = R(ivx,1) - bleftx(i+1,1)*2*etan(i,j)/xstp(j)/xstpc(j+1);

% Top Vx boundary condtion
L(c) = L(c) + btopx(j+1,2)*etas(i,j+1)/ystpc(i)/ystp(i);
R(ivx,1)=R(ivx,1)-btopx(j+1,1)*etas(i,j+1)/ystpc(i)/ystp(i);

% Right Vx nodes
c = c + 1;
I(c) = ivx;
J(c) = ivx + ynum3;
L(c) = 2*etan(i,j+1)/xstp(j+1)/xstpc(j+1);

% Bottom Vx nodes
c = c + 1;
I(c) = ivx;
J(c) = ivx + 3;
L(c) = etas(i+1,j+1)/ystpc(i+1)/ystp(i);

% Top Left Vy boundary condition applied to bottom left
c = c + 1;
I(c) = ivx;
J(c) = ivx + 1;
L(c) = btopy(j+1,2)*etas(i,j+1)/xstpc(j+1)/ystp(i);
R(ivx,1) = R(ivx,1) - btopy(j+1,1)*etas(i,j+1)/xstpc(j+1)/ystp(i);

% Bottom Left Vy node
L(c) = L(c) - etas(i+1,j+1)/xstpc(j+1)/ystp(i);

% Top Right Vy boundary condition applied to bottom right
c = c + 1;
I(c) = ivx;
J(c) = ivx + ynum3 + 1;
L(c) = -btopy(j+2,2)*etas(i,j+1)/xstpc(j+1)/ystp(i);
R(ivx,1) = R(ivx,1) + btopy(j+2,1)*etas(i,j+1)/xstpc(j+1)/ystp(i);

% Bottom Right Vy node
L(c) = L(c) + etas(i+1,j+1)/xstpc(j+1)/ystp(i);

% Left P node
c = c + 1;
I(c) = ivx;
J(c) = ivx + 2;
L(c) = pscale/xstpc(j+1);

% Right P node
c = c + 1;
I(c) = ivx;
J(c) = ivx + ynum3 + 2;
L(c) = -pscale/xstpc(j+1);

% y-Stokes coefficients

% Right Part
R(ivy,1) = RY(i+1,j+1);

% Central Vy node
c = c + 1;
I(c) = ivy;
J(c) = ivy;
L(c) = -2*(etan(i+1,j)/ystp(i+1)+etan(i,j)/ystp(i))/ystpc(i+1)-(etas(i+1,j+1)/xstpc(j+1)+etas(i+1,j)/xstpc(j))/xstp(j);

% Top Vy node
L(c) = L(c) + btopy(j+1,2)*2*etan(i,j)/ystp(i)/ystpc(i+1);
R(ivy,1) = R(ivy,1) - btopy(j+1,1)*2*etan(i,j)/ystp(i)/ystpc(i+1);

% Left Vy node add boundary condition for the left Vy node
L(c) = L(c) + blefty(i+1,2)*etas(i+1,j)/xstpc(j)/xstp(j);
R(ivy,1) = R(ivy,1) - blefty(i+1,1)*etas(i+1,j)/xstpc(j)/xstp(j);

% Bottom Vy node
c = c + 1;
I(c) = ivy;
J(c) = ivy + 3;
L(c) = 2*etan(i+1,j)/ystp(i+1)/ystpc(i+1);

% Right Vy node
c = c + 1;
I(c) = ivy;
J(c) = ivy + ynum3;
L(c) = etas(i+1,j+1)/xstpc(j+1)/xstp(j);

% Top left Vx node, boundary condition applied to top right Vx
c = c + 1;
I(c) = ivy;
J(c) = ivy - 1;
L(c) = bleftx(i+1,2)*etas(i+1,j)/ystpc(i+1)/xstp(j);
R(ivy,1) = R(ivy,1) - bleftx(i+1,1)*etas(i+1,j)/ystpc(i+1)/xstp(j);

% Top right Vx node using top Left Vx node
L(c) = L(c) - etas(i+1,j+1)/ystpc(i+1)/xstp(j);

% Bottom left Vx node, boundary condition applied to bottom right Vx
c = c + 1;
I(c) = ivy;
J(c) = ivy - 1 + 3;
L(c) = -bleftx(i+2,2)*etas(i+1,j)/ystpc(i+1)/xstp(j);
R(ivy,1) = R(ivy,1) + bleftx(i+2,1)*etas(i+1,j)/ystpc(i+1)/xstp(j);

% Bottom right Vx node using bottom Left Vx node
L(c) = L(c) + etas(i+1,j+1)/ystpc(i+1)/xstp(j);

% Top P node
c = c + 1;
I(c) = ivy;
J(c) = ivy + 1;
L(c) = pscale/ystpc(i+1);

% Bottom P node
c = c + 1;
I(c) = ivy;
J(c) = ivy + 1 + 3;
L(c) = -pscale/ystpc(i+1);

% Continuity equation

% Pressure definition for the boundary condition regions
c = c + 1;
I(c) = ipr;
J(c) = ipr;
% Pressure definition in one cell
L(c) = 2*pscale/(xstpavr+ystpavr);
R(ipr,1) = 2*prnorm/(xstpavr+ystpavr);

%------------------------------------------------------------------
%% Left boundary
%------------------------------------------------------------------

j = 1;
for i = 2:ynum-2
    
    % Indexes for P,vx,vy
    ivx=((j-1)*(ynum-1)+(i-1))*3+1;
    ivy=ivx+1;
    ipr=ivx+2;
    
    % x-Stokes coefficients
    
    % Right hand side
    R(ivx,1) = RX(i+1,j+1);
    
    % Central Vx node
    c = c + 1;
    I(c) = ivx;
    J(c) = ivx;
    L(c) = -2*(etan(i,j+1)/xstp(j+1)+etan(i,j)/xstp(j))/xstpc(j+1) ...
        -(etas(i+1,j+1)/ystpc(i+1)+etas(i,j+1)/ystpc(i))/ystp(i);
    
    % Left Vx boundary condition
    L(c) = L(c) + bleftx(i+1,2)*2*etan(i,j)/xstp(j)/xstpc(j+1);
    R(ivx,1) = R(ivx,1) - bleftx(i+1,1)*2*etan(i,j)/xstp(j)/xstpc(j+1);
    
    % Right Vx nodes
    c = c + 1;
    I(c) = ivx;
    J(c) = ivx + ynum3;
    L(c) = 2*etan(i,j+1)/xstp(j+1)/xstpc(j+1);
    
    % Top Vx node
    c = c + 1;
    I(c) = ivx;
    J(c) = ivx - 3;
    L(c) = etas(i,j+1)/ystpc(i)/ystp(i);
    
    % Bottom Vx nodes
    c = c + 1;
    I(c) = ivx;
    J(c) = ivx + 3;
    L(c) = etas(i+1,j+1)/ystpc(i+1)/ystp(i);
    
    % Top Left Vy node
    c = c + 1;
    I(c) = ivx;
    J(c) = ivx - 3 + 1;
    L(c) = etas(i,j+1)/xstpc(j+1)/ystp(i);
    
    % Top Right Vy node
    c = c + 1;
    I(c) = ivx;
    J(c) = ivx - 3 + 1 + ynum3;
    L(c) = -etas(i,j+1)/xstpc(j+1)/ystp(i);
    
    % Bottom Left Vy node
    c = c + 1;
    I(c) = ivx;
    J(c) = ivx + 1;
    L(c) = -etas(i+1,j+1)/xstpc(j+1)/ystp(i);
    
    % Bottom Right Vy node
    c = c + 1;
    I(c) = ivx;
    J(c) = ivx + 1 + ynum3;
    L(c) = etas(i+1,j+1)/xstpc(j+1)/ystp(i);
    
    % Left P node
    c = c + 1;
    I(c) = ivx;
    J(c) = ivx + 2;
    L(c) = pscale/xstpc(j+1);
    
    % Right P node
    c = c + 1;
    I(c) = ivx;
    J(c) = ivx + ynum3 + 2;
    L(c) = -pscale/xstpc(j+1);
    
    % y-Stokes coefficients
    
    % Right Part
    R(ivy,1) = RY(i+1,j+1);
    
    % Central Vy node
    c = c + 1;
    I(c) = ivy;
    J(c) = ivy;
    L(c) = -2*(etan(i+1,j)/ystp(i+1)+etan(i,j)/ystp(i))/ystpc(i+1)-(etas(i+1,j+1)/xstpc(j+1)+etas(i+1,j)/xstpc(j))/xstp(j);
    
    % Left Vy node add boundary condition for the left Vy node
    L(c) = L(c) + blefty(i+1,2)*etas(i+1,j)/xstpc(j)/xstp(j);
    R(ivy,1)= R(ivy,1) - blefty(i+1,1)*etas(i+1,j)/xstpc(j)/xstp(j);
  
    % Bottom Vy node
    if i < ynum-2
        c = c + 1;
        I(c) = ivy;
        J(c) = ivy + 3;
        L(c) = 2*etan(i+1,j)/ystp(i+1)/ystpc(i+1);
    else
        % Add boundary condition for the bottom Vy node
        L(c) = L(c) + bbottomy(j+1,2)*2*etan(i+1,j)/ystp(i+1)/ystpc(i+1);
        R(ivy,1) = R(ivy,1) - bbottomy(j+1,1)*2*etan(i+1,j)/ystp(i+1)/ystpc(i+1);
    end
 
    % Top Vy node
    c = c + 1;
    I(c) = ivy;
    J(c) = ivy - 3;
    L(c) = 2*etan(i,j)/ystp(i)/ystpc(i+1);
    
    % Right Vy node
    c = c + 1;
    I(c) = ivy;
    J(c) = ivy + ynum3;
    L(c) = etas(i+1,j+1)/xstpc(j+1)/xstp(j);
    
    % Top left Vx node, boundary condition applied to top right Vx
    c = c + 1;
    I(c) = ivy;
    J(c) = ivy - 1;
    L(c) = bleftx(i+1,2)*etas(i+1,j)/ystpc(i+1)/xstp(j);
    R(ivy,1) = R(ivy,1) - bleftx(i+1,1)*etas(i+1,j)/ystpc(i+1)/xstp(j);
    
    % Top right Vx node using top Left Vx node
    L(c) = L(c) - etas(i+1,j+1)/ystpc(i+1)/xstp(j);
    
    % Bottom left Vx node, boundary condition applied to bottom right Vx
    c = c + 1;
    I(c) = ivy;
    J(c) = ivy - 1 + 3;
    L(c) = -bleftx(i+2,2)*etas(i+1,j)/ystpc(i+1)/xstp(j);
    R(ivy,1) = R(ivy,1) + bleftx(i+2,1)*etas(i+1,j)/ystpc(i+1)/xstp(j);
    
    % Bottom right Vx node using bottom Left Vx node
    L(c) = L(c) + etas(i+1,j+1)/ystpc(i+1)/xstp(j);
    
    % Top P node
    c = c + 1;
    I(c) = ivy;
    J(c) = ivy + 1;
    L(c) = pscale/ystpc(i+1);
    
    % Bottom P node
    c = c + 1;
    I(c) = ivy;
    J(c) = ivy + 1 + 3;
    L(c) = -pscale/ystpc(i+1);
    
    % Continuity equation
    
    % Right Part
    R(ipr,1)=RC(i,j);
    
    % Computing Current Continuity coefficients
    
    % Right Vx node
    c = c + 1;
    I(c) = ipr;
    J(c) = ipr - 2;
    L(c) = pscale/xstp(j);
    
    % Add boundary condition for the left Vx node
    L(c) = L(c) - bleftx(i+1,2)*pscale/xstp(j);
    R(ipr,1) = R(ipr,1) + bleftx(i+1,1)*pscale/xstp(j);

    % Top Vy node
    c = c + 1;
    I(c) = ipr;
    J(c) = ipr -1 - 3;
    L(c) = -pscale/ystp(i);

    % Bottom Vy node
    c = c + 1;
    I(c) = ipr;
    J(c) = ipr - 1;
    L(c) = pscale/ystp(i);
     
end
    
%------------------------------------------------------------------
%% Bottom left corner
%------------------------------------------------------------------
i = ynum-1;
j = 1;
% Indexes for P,vx,vy
ivx = ((j-1)*(ynum-1)+(i-1))*3+1;
ivy = ivx+1;
ipr = ivx+2;

% x-Stokes coefficients
R(ivx,1) = RX(i+1,j+1);

% Central Vx node
c = c + 1;
I(c) = ivx;
J(c) = ivx;
L(c) = -2*(etan(i,j+1)/xstp(j+1)+etan(i,j)/xstp(j))/xstpc(j+1)-(etas(i+1,j+1)/ystpc(i+1)+etas(i,j+1)/ystpc(i))/ystp(i);
 
% Left Vx boundary condition
L(c) = L(c) + bleftx(i+1,2)*2*etan(i,j)/xstp(j)/xstpc(j+1);
R(ivx,1) = R(ivx,1)-bleftx(i+1,1)*2*etan(i,j)/xstp(j)/xstpc(j+1);

% Bottom Vx boundary condition
L(c) = L(c) + bbottomx(j+1,2)*etas(i+1,j+1)/ystpc(i+1)/ystp(i);
R(ivx,1) = R(ivx,1) - bbottomx(j+1,1)*etas(i+1,j+1)/ystpc(i+1)/ystp(i);

% Right Vx node
c = c + 1;
I(c) = ivx;
J(c) = ivx + ynum3;
L(c) = 2*etan(i,j+1)/xstp(j+1)/xstpc(j+1);

% Top Vx node
c = c + 1;
I(c) = ivx;
J(c) = ivx - 3;
L(c) = etas(i,j+1)/ystpc(i)/ystp(i);

% Top left Vy node
c = c + 1;
I(c) = ivx;
J(c) = ivx - 3 + 1;
L(c) = etas(i,j+1)/xstpc(j+1)/ystp(i);

% Bottom left Vy node, apply boundary condition
L(c) = L(c) - bbottomy(j+1,2)*etas(i+1,j+1)/xstpc(j+1)/ystp(i);
R(ivx,1) = R(ivx,1) + bbottomy(j+1,1)*etas(i+1,j+1)/xstpc(j+1)/ystp(i);
           
% Top right Vy node
c = c + 1;
I(c) = ivx;
J(c) = ivx - 3 + 1 + ynum3;
L(c) = -etas(i,j+1)/xstpc(j+1)/ystp(i);

% Bottom right Vy node, apply boundary condition
L(c) = L(c) + bbottomy(j+2,2)*etas(i+1,j+1)/xstpc(j+1)/ystp(i);
R(ivx,1) = R(ivx,1) - bbottomy(j+2,1)*etas(i+1,j+1)/xstpc(j+1)/ystp(i);
 
% Left P node
c = c + 1;
I(c) = ivx;
J(c) = ivx + 2;
L(c) = pscale/xstpc(j+1);

% Right P node
c = c + 1;
I(c) = ivx;
J(c) = ivx + 2 + ynum3;
L(c) = -pscale/xstpc(j+1);

% y- Stokes coefficients
% Ghost Vy-parameter=0 for numbering  (i = ynum-1)

% Central Vy node
c = c + 1;
I(c) = ivy;
J(c) = ivy;
L(c) = 2*pscale/(xstpavr+ystpavr);

% Continuity equation

% Right Part
R(ipr,1) = RC(i,j);

% Right Vx node
c = c + 1;
I(c) = ipr;
J(c) = ipr - 2 ;
L(c) = pscale/xstp(j);

% Left Vx node, apply boundary condition
L(c) = L(c) - bleftx(i+1,2)*pscale/xstp(j);
R(ipr,1) = R(ipr,1) + bleftx(i+1,1)*pscale/xstp(j);

% Top Vy node
c = c + 1;
I(c) = ipr;
J(c) = ipr - 1 - 3;
L(c)=-pscale/ystp(i);

% Bottom Vy node, apply boundary condition
L(c) = L(c) + bbottomy(j+1,2)*pscale/ystp(i);
R(ipr,1) = R(ipr,1) - bbottomy(j+1,1)*pscale/ystp(i);

%------------------------------------------------------------------
%% Top boundary
%------------------------------------------------------------------
i = 1;

for j = 2:xnum-2
    
    % Indexes for P,vx,vy
    ivx=((j-1)*(ynum-1)+(i-1))*3+1;
    ivy=ivx+1;
    ipr=ivx+2;
    
    % x-Stokes coefficients
    
    % Right hand side
    R(ivx,1) = RX(i+1,j+1);
    
    % Central Vx node
    c = c + 1;
    I(c) = ivx;
    J(c) = ivx;
    L(c) = -2*(etan(i,j+1)/xstp(j+1)+etan(i,j)/xstp(j))/xstpc(j+1) ...
        -(etas(i+1,j+1)/ystpc(i+1)+etas(i,j+1)/ystpc(i))/ystp(i);
    
    % Top Vx, apply boundary conition
    L(c) = L(c) + btopx(j+1,2)*etas(i,j+1)/ystpc(i)/ystp(i);
    R(ivx,1) = R(ivx,1) - btopx(j+1,1)*etas(i,j+1)/ystpc(i)/ystp(i);
    
    % Right Vx node
    if j < xnum-2
        c = c + 1;
        I(c) = ivx;
        J(c) = ivx + ynum3;
        L(c) = 2*etan(i,j+1)/xstp(j+1)/xstpc(j+1);
    else
        L(c) = L(c) + brightx(i+1,2)*2*etan(i,j+1)/xstp(j+1)/xstpc(j+1);
        R(ivx,1) = R(ivx,1) - brightx(i+1,1)*2*etan(i,j+1)/xstp(j+1)/xstpc(j+1);
    end
    
    % Left Vx node
    c = c + 1;
    I(c) = ivx;
    J(c) = ivx - ynum3;
    L(c) = 2*etan(i,j)/xstp(j)/xstpc(j+1);

    % Bottom Vx nodes
    c = c + 1;
    I(c) = ivx;
    J(c) = ivx + 3;
    L(c) = etas(i+1,j+1)/ystpc(i+1)/ystp(i);
    
    % Top Left Vy node, apply boundary condition to bottom left Vy node
    c = c + 1;
    I(c) = ivx;
    J(c) = ivx + 1;
    L(c) = btopy(j+1,2)*etas(i,j+1)/xstpc(j+1)/ystp(i);
    R(ivx,1) = R(ivx,1) - btopy(j+1,1)*etas(i,j+1)/xstpc(j+1)/ystp(i);
    
    % Bottom left Vy node, apply boundary condition
    L(c) = L(c) - etas(i+1,j+1)/xstpc(j+1)/ystp(i);

    % Top Right Vy node, apply boundary condttion to bottom right Vy node
    c = c + 1;
    I(c) = ivx;
    J(c) = ivx + 1 + ynum3;
    L(c) = -btopy(j+2,2)*etas(i,j+1)/xstpc(j+1)/ystp(i);
    R(ivx,1) = R(ivx,1) + btopy(j+2,1)*etas(i,j+1)/xstpc(j+1)/ystp(i);
    
    % Bottom Right Vy node, apply boundary condition
    L(c) = L(c) + etas(i+1,j+1)/xstpc(j+1)/ystp(i);
    
    % Left P node
    c = c + 1;
    I(c) = ivx;
    J(c) = ivx + 2;
    L(c) = pscale/xstpc(j+1);
    
    % Right P node
    c = c + 1;
    I(c) = ivx;
    J(c) = ivx + ynum3 + 2;
    L(c) = -pscale/xstpc(j+1);
    
    % y-Stokes coefficients
    
    % Right Part
    R(ivy,1) = RY(i+1,j+1);
    
    % Central Vy node
    c = c + 1;
    I(c) = ivy;
    J(c) = ivy;
    L(c) = -2*(etan(i+1,j)/ystp(i+1)+etan(i,j)/ystp(i))/ystpc(i+1)-(etas(i+1,j+1)/xstpc(j+1)+etas(i+1,j)/xstpc(j))/xstp(j);
    
    % Add boundary condition for the top Vy node
    L(c) = L(c) + btopy(j+1,2)*2*etan(i,j)/ystp(i)/ystpc(i+1);
    R(ivy,1) = R(ivy,1) - btopy(j+1,1)*2*etan(i,j)/ystp(i)/ystpc(i+1);
        
    % Bottom Vy node
    c = c + 1;
    I(c) = ivy;
    J(c) = ivy + 3;
    L(c) = 2*etan(i+1,j)/ystp(i+1)/ystpc(i+1);
    
    % Left Vy node
    c = c + 1;
    I(c) = ivy;
    J(c) = ivy - ynum3;
    L(c) = etas(i+1,j)/xstpc(j)/xstp(j);
    
    % Right Vy node
    c = c + 1;
    I(c) = ivy;
    J(c) = ivy + ynum3;
    L(c) = etas(i+1,j+1)/xstpc(j+1)/xstp(j);
    
    % Top left Vx node
    c = c + 1;
    I(c) = ivy;
    J(c) = ivy - 1 - ynum3;
    L(c) = etas(i+1,j)/ystpc(i+1)/xstp(j);
    
    % Top right Vx node
    c = c + 1;
    I(c) = ivy;
    J(c) = ivy - 1;
    L(c)=-etas(i+1,j+1)/ystpc(i+1)/xstp(j);
    
    % Bottom left Vx node
    c = c + 1;
    I(c) = ivy;
    J(c) = ivy - 1 + 3 - ynum3;
    L(c) = -etas(i+1,j)/ystpc(i+1)/xstp(j);
     
    % Bottom right Vx node
    c = c + 1;
    I(c) = ivy;
    J(c) = ivy - 1 + 3;
    L(c) = etas(i+1,j+1)/ystpc(i+1)/xstp(j);
    
    % Top P node
    c = c + 1;
    I(c) = ivy;
    J(c) = ivy + 1;
    L(c) = pscale/ystpc(i+1);
    
    % Bottom P node
    c = c + 1;
    I(c) = ivy;
    J(c) = ivy + 1 + 3;
    L(c) = -pscale/ystpc(i+1);
    
    % Continuity equation
    
    % Right Part
    R(ipr,1)=RC(i,j);
    
    % Computing Current Continuity coefficients
    
    % Left Vx node
    c = c + 1;
    I(c) = ipr;
    J(c) = ipr -2 - ynum3;
    L(c) = -pscale/xstp(j);
    
    % Right Vx node
    c = c + 1;
    I(c) = ipr;
    J(c) = ipr - 2;
    L(c) = pscale/xstp(j);

    % Bottom Vy node
    c = c + 1;
    I(c) = ipr;
    J(c) = ipr - 1;
    L(c) = pscale/ystp(i);
    
    % Add boundary condition for the top Vy node
    L(c) = L(c) - btopy(j+1,2)*pscale/ystp(i);
    R(ipr,1) = R(ipr,1) + btopy(j+1,1)*pscale/ystp(i);
end


%------------------------------------------------------------------
%% Internal boundary conditions
%------------------------------------------------------------------

% x-Stokes coefficients internal boundary conditions
if bintern(1) > 0
    
    j = bintern(1);
    i = bintern(2):bintern(3);
    ivx_bounds = ((j-1)*(ynum-1)+(i-1))*3+1;
    
    for k = ivx_bounds
        c = c + 1;
        I(c) = k;
        J(c) = k;
        L(c) = 2*pscale/(xstpavr+ystpavr);
        R(k) = 2*pscale/(xstpavr+ystpavr)*bintern(4);
    end
    
else
    ivx_bounds = -1;
    
end

% y-Stokes coefficients internal boundary conditions
if bintern(5) > 0
    
    i = bintern(5);
    j = bintern(6):bintern(7);
    ivx = ((j-1)*(ynum-1)+(i-1))*3+1;
    ivy_bounds = ivx+1;
    for k = ivy_bounds
        c = c + 1;
        I(c) = k;
        J(c) = k;
        L(c) = 2*pscale/(xstpavr+ystpavr);
        R(k) = 2*pscale/(xstpavr+ystpavr)*bintern(8);
    end
    
else    
    ivy_bounds = -1;
    
end

%------------------------------------------------------------------
%% Internal points, no boundary conditions to apply, except
%  ghost nodes at i = ynum-2 and j = xnum-2
%------------------------------------------------------------------

for i = 2:ynum-2
    for j = 2:xnum-2
            
        % Indexes for P,vx,vy
        ivx=((j-1)*(ynum-1)+(i-1))*3+1;
        ivy=ivx+1;
        ipr=ivx+2;
        
        % x-Stokes coefficients
        if ~any(ivx == ivx_bounds) % Exclude internal boundary nodes
            
            % Central Vx node
            c = c + 1;
            I(c) = ivx;
            J(c) = ivx;
            
            % Right hand side
            R(ivx,1) = RX(i+1,j+1);
            
            L(c) = -2*(etan(i,j+1)/xstp(j+1)+etan(i,j)/xstp(j))/xstpc(j+1) ...
                -(etas(i+1,j+1)/ystpc(i+1)+etas(i,j+1)/ystpc(i))/ystp(i);
            
            % Right Vx node
            if j < xnum-2
                c = c + 1;
                I(c) = ivx;
                J(c) = ivx + ynum3;
                L(c) = 2*etan(i,j+1)/xstp(j+1)/xstpc(j+1);
            else
                L(c) = L(c) + brightx(i+1,2)*2*etan(i,j+1)/xstp(j+1)/xstpc(j+1);
                R(ivx,1) = R(ivx,1) - brightx(i+1,1)*2*etan(i,j+1)/xstp(j+1)/xstpc(j+1);
            end
            
            % Left Vx node
            c = c + 1;
            I(c) = ivx;
            J(c) = ivx - ynum3;
            L(c) = 2*etan(i,j)/xstp(j)/xstpc(j+1);
            
            % Top Vx node
            c = c + 1;
            I(c) = ivx;
            J(c) = ivx - 3;
            L(c) = etas(i,j+1)/ystpc(i)/ystp(i);
            
            % Bottom Vx node
            c = c + 1;
            I(c) = ivx;
            J(c) = ivx + 3;
            L(c) = etas(i+1,j+1)/ystpc(i+1)/ystp(i);
            
            % Top Left Vy node
            c = c + 1;
            I(c) = ivx;
            J(c) = ivx - 3 + 1;
            L(c) = etas(i,j+1)/xstpc(j+1)/ystp(i);
            
            % Top right Vy node
            c = c + 1;
            I(c) = ivx;
            J(c) = ivx - 3 + 1 + ynum3;
            L(c) = -etas(i,j+1)/xstpc(j+1)/ystp(i);
            
            % Bottom left Vy node
            c = c + 1;
            I(c) = ivx;
            J(c) = ivx + 1;
            L(c) = -etas(i+1,j+1)/xstpc(j+1)/ystp(i);
            
            % Bottom Right Vy node
            c = c + 1;
            I(c) = ivx;
            J(c) = ivx + 1 + ynum3;
            L(c) = etas(i+1,j+1)/xstpc(j+1)/ystp(i);
            
            % Left P node
            c = c + 1;
            I(c) = ivx;
            J(c) = ivx + 2;
            L(c) = pscale/xstpc(j+1);
            
            % Right P node
            c = c + 1;
            I(c) = ivx;
            J(c) = ivx + ynum3 + 2;
            L(c) = -pscale/xstpc(j+1);
        end
        
        % y-Stokes coefficients
        if ~any(ivy == ivy_bounds) % Exclude internal boundary nodes
            
            % Right Part
            R(ivy,1) = RY(i+1,j+1);
            
            % Central Vy node
            c = c + 1;
            I(c) = ivy;
            J(c) = ivy;
            L(c) = -2*(etan(i+1,j)/ystp(i+1)+etan(i,j)/ystp(i))/ystpc(i+1)-(etas(i+1,j+1)/xstpc(j+1)+etas(i+1,j)/xstpc(j))/xstp(j);
            
            % Bottom Vy node
            if i < ynum-2
                c = c + 1;
                I(c) = ivy;
                J(c) = ivy + 3;
                L(c) = 2*etan(i+1,j)/ystp(i+1)/ystpc(i+1);
            else
                L(c) = L(c) + bbottomy(j+1,2)*2*etan(i+1,j)/ystp(i+1)/ystpc(i+1);
                R(ivy,1) = R(ivy,1) - bbottomy(j+1,1)*2*etan(i+1,j)/ystp(i+1)/ystpc(i+1);
            end
            
            % Top Vy node
            c = c + 1;
            I(c) = ivy;
            J(c) = ivy - 3;
            L(c) = 2*etan(i,j)/ystp(i)/ystpc(i+1);
            
            % Left Vy node
            c = c + 1;
            I(c) = ivy;
            J(c) = ivy - ynum3;
            L(c) = etas(i+1,j)/xstpc(j)/xstp(j);
            
            % Right Vy node
            c = c + 1;
            I(c) = ivy;
            J(c) = ivy + ynum3;
            L(c) = etas(i+1,j+1)/xstpc(j+1)/xstp(j);
            
            % Top left Vx node
            c = c + 1;
            I(c) = ivy;
            J(c) = ivy - 1 - ynum3;
            L(c) = etas(i+1,j)/ystpc(i+1)/xstp(j);
            
            % Top right Vx node
            c = c + 1;
            I(c) = ivy;
            J(c) = ivy - 1;
            L(c) = -etas(i+1,j+1)/ystpc(i+1)/xstp(j);
            
            % Bottom left Vx node
            c = c + 1;
            I(c) = ivy;
            J(c) = ivy - 1 + 3 - ynum3;
            L(c) = -etas(i+1,j)/ystpc(i+1)/xstp(j);
            
            % Bottom right Vx node
            c = c + 1;
            I(c) = ivy;
            J(c) = ivy - 1 + 3;
            L(c) = etas(i+1,j+1)/ystpc(i+1)/xstp(j);
            
            % Top P node
            c = c + 1;
            I(c) = ivy;
            J(c) = ivy + 1;
            L(c) = pscale/ystpc(i+1);
            
            % Bottom P node
            c = c + 1;
            I(c) = ivy;
            J(c) = ivy + 1 + 3;
            L(c) = -pscale/ystpc(i+1);
            
            % Continuity equation
            
            % Right Part
            R(ipr,1)=RC(i,j);
            
            % Computing Current Continuity coefficients
            
            % Left Vx node
            c = c + 1;
            I(c) = ipr;
            J(c) = ipr -2 - ynum3;
            L(c) = -pscale/xstp(j);
            
            % Right Vx node
            c = c + 1;
            I(c) = ipr;
            J(c) = ipr - 2;
            L(c) = pscale/xstp(j);
            
            % Top Vy node
            c = c + 1;
            I(c) = ipr;
            J(c) = ipr - 1 - 3;
            L(c) = -pscale/ystp(i);
            
            % Bottom Vy node
            c = c + 1;
            I(c) = ipr;
            J(c) = ipr - 1;
            L(c) = pscale/ystp(i);
        end
        
    end
end

%------------------------------------------------------------------
%% Bottom boundary
%------------------------------------------------------------------
i = ynum-1;

for j = 2:xnum-2
    
    % Indexes for P,vx,vy
    ivx=((j-1)*(ynum-1)+(i-1))*3+1;
    ivy=ivx+1;
    ipr=ivx+2;
    
    % x-Stokes coefficients
    
    % Right hand side
    R(ivx,1) = RX(i+1,j+1);
    
    % Central Vx node
    c = c + 1;
    I(c) = ivx;
    J(c) = ivx;
    L(c) = -2*(etan(i,j+1)/xstp(j+1)+etan(i,j)/xstp(j))/xstpc(j+1) ...
        -(etas(i+1,j+1)/ystpc(i+1)+etas(i,j+1)/ystpc(i))/ystp(i);
    
    % Bottom Vx boundary condition
    L(c) = L(c) + bbottomx(j+1,2)*etas(i+1,j+1)/ystpc(i+1)/ystp(i);
    R(ivx,1) = R(ivx,1) - bbottomx(j+1,1)*etas(i+1,j+1)/ystpc(i+1)/ystp(i);

    % Right Vx node, ghost nodes BC
    if j< xnum-2
        c = c + 1;
        I(c) = ivx;
        J(c) = ivx + ynum3;
        L(c) = 2*etan(i,j+1)/xstp(j+1)/xstpc(j+1);
    else
        L(c) = L(c) + brightx(i+1,2)*2*etan(i,j+1)/xstp(j+1)/xstpc(j+1);
        R(ivx,1) = R(ivx,1) - brightx(i+1,1)*2*etan(i,j+1)/xstp(j+1)/xstpc(j+1);
    end
    
    % Left Vx node
    c = c + 1;
    I(c) = ivx;
    J(c) = ivx - ynum3;
    L(c) = 2*etan(i,j)/xstp(j)/xstpc(j+1);
     
    % Top Vx node
    c = c + 1;
    I(c) = ivx;
    J(c) = ivx - 3;
    L(c) = etas(i,j+1)/ystpc(i)/ystp(i);
    
    % Top left Vy node
    c = c + 1;
    I(c) = ivx;
    J(c) = ivx - 3 + 1;
    L(c) = etas(i,j+1)/xstpc(j+1)/ystp(i);
    
    % Bottom left Vy node, apply boundary condition
    L(c) = L(c) - bbottomy(j+1,2)*etas(i+1,j+1)/xstpc(j+1)/ystp(i);
    R(ivx,1) = R(ivx,1) + bbottomy(j+1,1)*etas(i+1,j+1)/xstpc(j+1)/ystp(i);
    
    % Top right Vy node
    c = c + 1;
    I(c) = ivx;
    J(c) = ivx - 3 + 1 + ynum3;
    L(c) = -etas(i,j+1)/xstpc(j+1)/ystp(i);
    
    % Bottom right Vy node, apply boundary condition
    L(c) = L(c) + bbottomy(j+2,2)*etas(i+1,j+1)/xstpc(j+1)/ystp(i);
    R(ivx,1) = R(ivx,1) - bbottomy(j+2,1)*etas(i+1,j+1)/xstpc(j+1)/ystp(i);
    
    % Left P node
    c = c + 1;
    I(c) = ivx;
    J(c) = ivx + 2;
    L(c) = pscale/xstpc(j+1);
    
    % Right P node
    c = c + 1;
    I(c) = ivx;
    J(c) = ivx + 2 + ynum3;
    L(c) = -pscale/xstpc(j+1);
    
    % y- Stokes coefficients
    % Ghost Vy_parameter = 0 for numbering
    
    % Central Vy node
    c = c + 1;
    I(c) = ivy;
    J(c) = ivy;
    L(c) = 2*pscale/(xstpavr+ystpavr);
    R(ivy,1) = 0;
    
    % Continuity equation
    
    % Right Part
    R(ipr,1) = RC(i,j);
    
    % Left Vx node
    c = c + 1;
    I(c) = ipr;
    J(c) = ipr - 2 - ynum3;
    L(c) = -pscale/xstp(j);
    
    % Right Vx node
    c = c + 1;
    I(c) = ipr;
    J(c) = ipr - 2 ;
    L(c) = pscale/xstp(j);
    
    % Top Vy node
    c = c + 1;
    I(c) = ipr;
    J(c) = ipr - 1 - 3;
    L(c) = -pscale/ystp(i);
    
    % Bottom Vy node, apply boundary condition
    L(c) = L(c) + bbottomy(j+1,2)*pscale/ystp(i);
    R(ipr,1) = R(ipr,1) - bbottomy(j+1,1)*pscale/ystp(i); 
    
end

%------------------------------------------------------------------
%% Top right corner
%------------------------------------------------------------------
i = 1;
j = xnum -1;

% Indexes for P,vx,vy
ivx=((j-1)*(ynum-1)+(i-1))*3+1;
ivy=ivx+1;
ipr=ivx+2;

% x-Stokes coefficients
% Ghost Vx_parameter = 0 for numbering

% Central Vx node
c = c + 1;
I(c) = ivx;
J(c) = ivx;
L(c) = 2*pscale/(xstpavr+ystpavr);
R(ivx,1) = 0;

% y-Stokes coefficients

% Right Part
R(ivy,1) = RY(i+1,j+1);

% Central Vy node
c = c + 1;
I(c) = ivy;
J(c) = ivy;
L(c) = -2*(etan(i+1,j)/ystp(i+1)+etan(i,j)/ystp(i))/ystpc(i+1)-(etas(i+1,j+1)/xstpc(j+1)+etas(i+1,j)/xstpc(j))/xstp(j);

% Add boundary condition for the top Vy node
L(c) = L(c) + btopy(j+1,2)*2*etan(i,j)/ystp(i)/ystpc(i+1);
R(ivy,1) = R(ivy,1) - btopy(j+1,1)*2*etan(i,j)/ystp(i)/ystpc(i+1);

% Right Vy node, add boundary condition
L(c) = L(c) + brighty(i+1,2)*etas(i+1,j+1)/xstpc(j+1)/xstp(j);
R(ivy,1) = R(ivy,1) - brighty(i+1,1)*etas(i+1,j+1)/xstpc(j+1)/xstp(j);

% Bottom Vy node
c = c + 1;
I(c) = ivy;
J(c) = ivy + 3;
L(c) = 2*etan(i+1,j)/ystp(i+1)/ystpc(i+1);

% Left Vy node
c = c + 1;
I(c) = ivy;
J(c) = ivy - ynum3;
L(c) = etas(i+1,j)/xstpc(j)/xstp(j);

% Top left Vx node
c = c + 1;
I(c) = ivy;
J(c) = ivy - 1 - ynum3;
L(c) = etas(i+1,j)/ystpc(i+1)/xstp(j);

% Top right Vx node, apply boundary condition to top left Vx
L(c) = L(c) - brightx(i+1,2)*etas(i+1,j+1)/ystpc(i+1)/xstp(j);
R(ivy,1) = R(ivy,1) + brightx(i+1,1)*etas(i+1,j+1)/ystpc(i+1)/xstp(j);

% Bottom left Vx node
c = c + 1;
I(c) = ivy;
J(c) = ivy - 1 + 3 - ynum3;
L(c) = -etas(i+1,j)/ystpc(i+1)/xstp(j);

% Bottom right Vx node, apply boundary condition to bottom left Vx
L(c) = L(c) + brightx(i+2,2)*etas(i+1,j+1)/ystpc(i+1)/xstp(j);
R(ivy,1) = R(ivy,1) - brightx(i+2,1)*etas(i+1,j+1)/ystpc(i+1)/xstp(j);

% Top P node
c = c + 1;
I(c) = ivy;
J(c) = ivy + 1;
L(c) = pscale/ystpc(i+1);

% Bottom P node
c = c + 1;
I(c) = ivy;
J(c) = ivy + 1 + 3;
L(c) = -pscale/ystpc(i+1);

% Continuity equation

% Right Part
R(ipr,1)=RC(i,j);

% Computing Current Continuity coefficients

% Left Vx node
c = c + 1;
I(c) = ipr;
J(c) = ipr -2 - ynum3;
L(c) = -pscale/xstp(j);

% Right Vx node, apply boundary condition to left Vx node
L(c) = L(c)+brightx(i+1,2)*pscale/xstp(j);
R(ipr,1) = R(ipr,1) - brightx(i+1,1)*pscale/xstp(j);

% Bottom Vy node
c = c + 1;
I(c) = ipr;
J(c) = ipr - 1;
L(c) = pscale/ystp(i);

% Top Vy node, add boundary condition to bottom Vy node
L(c) = L(c) - btopy(j+1,2)*pscale/ystp(i);
R(ipr,1) = R(ipr,1) + btopy(j+1,1)*pscale/ystp(i);

%------------------------------------------------------------------
%% Right boundary
%------------------------------------------------------------------
j = xnum-1;
for i = 2:ynum-2
    
    % Indexes for P,vx,vy
    ivx=((j-1)*(ynum-1)+(i-1))*3+1;
    ivy=ivx+1;
    ipr=ivx+2;
    
    % x-Stokes coefficients
    % Ghost Vx_parameter = 0 for numbering
    
    % Central Vx node
    c = c + 1;
    I(c) = ivx;
    J(c) = ivx;
    L(c) = 2*pscale/(xstpavr+ystpavr);
    R(ivx,1) = 0;
    
    % y-Stokes coefficients
    
    % Right Part
    R(ivy,1) = RY(i+1,j+1);
    
    % Central Vy node
    c = c + 1;
    I(c) = ivy;
    J(c) = ivy;
    L(c) = -2*(etan(i+1,j)/ystp(i+1)+etan(i,j)/ystp(i))/ystpc(i+1)-(etas(i+1,j+1)/xstpc(j+1)+etas(i+1,j)/xstpc(j))/xstp(j);
    
    % Right Vy node, add boundary condition
    L(c) = L(c) + brighty(i+1,2)*etas(i+1,j+1)/xstpc(j+1)/xstp(j);
    R(ivy,1) = R(ivy,1) - brighty(i+1,1)*etas(i+1,j+1)/xstpc(j+1)/xstp(j);
    
    % Bottom Vy node
    if (i<ynum-2)
        c = c + 1;
        I(c) = ivy;
        J(c) = ivy + 3;
        L(c) = 2*etan(i+1,j)/ystp(i+1)/ystpc(i+1);
    else
        % Add ghost boundary condition for the bottom Vy node
        L(c) = L(c) + bbottomy(j+1,2)*2*etan(i+1,j)/ystp(i+1)/ystpc(i+1);
        R(ivy,1) = R(ivy,1) - bbottomy(j+1,1)*2*etan(i+1,j)/ystp(i+1)/ystpc(i+1);
    end
    
    % Top Vy node
    c = c + 1;
    I(c) = ivy;
    J(c) = ivy - 3;
    L(c)=2*etan(i,j)/ystp(i)/ystpc(i+1);
   
    % Left Vy node
    c = c + 1;
    I(c) = ivy;
    J(c) = ivy - ynum3;
    L(c) = etas(i+1,j)/xstpc(j)/xstp(j);
    
    % Top left Vx node
    c = c + 1;
    I(c) = ivy;
    J(c) = ivy - 1 - ynum3;
    L(c) = etas(i+1,j)/ystpc(i+1)/xstp(j);
    
    % Top right Vx node, apply boundary condition to top left Vx
    L(c) = L(c) - brightx(i+1,2)*etas(i+1,j+1)/ystpc(i+1)/xstp(j);
    R(ivy,1) = R(ivy,1) + brightx(i+1,1)*etas(i+1,j+1)/ystpc(i+1)/xstp(j);
    
    % Bottom left Vx node
    c = c + 1;
    I(c) = ivy;
    J(c) = ivy - 1 + 3 - ynum3;
    L(c) = -etas(i+1,j)/ystpc(i+1)/xstp(j);
    
    % Bottom right Vx node, apply boundary condition to bottom left Vx
    L(c) = L(c) + brightx(i+2,2)*etas(i+1,j+1)/ystpc(i+1)/xstp(j);
    R(ivy,1) = R(ivy,1) - brightx(i+2,1)*etas(i+1,j+1)/ystpc(i+1)/xstp(j);
    
    % Top P node
    c = c + 1;
    I(c) = ivy;
    J(c) = ivy + 1;
    L(c) = pscale/ystpc(i+1);
    
    % Bottom P node
    c = c + 1;
    I(c) = ivy;
    J(c) = ivy + 1 + 3;
    L(c) = -pscale/ystpc(i+1);
    
    % Continuity equation
    
    % Right Part
    R(ipr,1)=RC(i,j);
    
    % Computing Current Continuity coefficients
    
    % Left Vx node
    c = c + 1;
    I(c) = ipr;
    J(c) = ipr -2 - ynum3;
    L(c) = -pscale/xstp(j);
    
    % Right Vx node, apply boundary condition to left Vx node
    L(c) = L(c)+brightx(i+1,2)*pscale/xstp(j);
    R(ipr,1) = R(ipr,1) - brightx(i+1,1)*pscale/xstp(j);
    
    % Top Vy node
    c = c + 1;
    I(c) = ipr;
    J(c) = ipr-1-3;
    L(c) = -pscale/ystp(i);
    
    % Bottom Vy node
    c = c + 1;
    I(c) = ipr;
    J(c) = ipr -1 ;
    L(c) = pscale/ystp(i);
    
end

%------------------------------------------------------------------
%% Bottom right corner
%------------------------------------------------------------------
i = ynum-1;
j = xnum-1;

% Indexes for P,vx,vy
ivx=((j-1)*(ynum-1)+(i-1))*3+1;
ivy=ivx+1;
ipr=ivx+2;

% x-Stokes coefficients
% Ghost Vx_parameter = 0 for numbering

% Central Vx node
c = c + 1;
I(c) = ivx;
J(c) = ivx;
L(c) = 2*pscale/(xstpavr+ystpavr);
R(ivx,1) = 0;

% y- Stokes coefficients
% Ghost Vy_parameter = 0 for numbering

% Central Vy node
c = c + 1;
I(c) = ivy;
J(c) = ivy;
L(c) = 2*pscale/(xstpavr+ystpavr);
R(ivy,1) = 0;

% Continuity equation

% Right Part
R(ipr,1)=RC(i,j);

% Computing Current Continuity coefficients

% Left Vx node
c = c + 1;
I(c) = ipr;
J(c) = ipr -2 - ynum3;
L(c) = -pscale/xstp(j);

% Right Vx node, apply boundary condition to left Vx node
L(c) = L(c)+brightx(i+1,2)*pscale/xstp(j);
R(ipr,1) = R(ipr,1) - brightx(i+1,1)*pscale/xstp(j);

% Top Vy node
c = c + 1;
I(c) = ipr;
J(c) = ipr-1-3;
L(c) = -pscale/ystp(i);

% Add boundary condition for the bottom Vy node
L(c) = L(c) + bbottomy(j+1,2)*pscale/ystp(i);
R(ipr,1) = R(ipr,1) - bbottomy(j+1,1)*pscale/ystp(i);


%------------------------------------------------------------------
%% Return
%------------------------------------------------------------------

I = I(1:c);
J = J(1:c);
L = L(1:c);

