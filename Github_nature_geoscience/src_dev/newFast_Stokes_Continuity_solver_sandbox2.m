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
function[vx,resx,vy,resy,pr,resc,resg]=newFast_Stokes_Continuity_solver_sandbox2(prfirst,etas,etan,xnum,ynum,gridx,gridy,RX,RY,RC,bleft,bright,btop,bbottom,bintern,useSuiteSparse)
%
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
% are used for boundary conditions see Figure 7.17 in Gerya's Book
%
% Koefficient for scaling pressure
% Gerya's pressure scaling:
% Average x and y steps

 xstpavr=(gridx(xnum)-gridx(1))/(xnum-1);
 ystpavr=(gridy(ynum)-gridy(1))/(ynum-1);
 pscale=2*etan(1)/(xstpavr+ystpavr);
% The following scaling was found to yield the lowest L2 matrix condition
% number, but higher residual???
% pscale = abs(2*max(etan(:))/(min(diff(gridx)) + min(diff(gridy))));

%tic;
%[I,J,L,R]=mexFastStokesContinuityMatrix2(prfirst,etas,etan,gridx,gridy,RX,RY,RC,bleft,bright,btop,bbottom,bintern,pscale);
%toc;

[I,J,L,R]=newFastStokesContinuityMatrix2(prfirst,etas,etan,xnum,ynum,gridx,gridy,RX,RY,RC,bleft,bright,btop,bbottom,bintern,pscale);

% Assamble sparse matrix and solve

if useSuiteSparse
    % SuiteSparse package sparse equivalent
    L2 = sparse2(I,J,L,(xnum-1)*(ynum-1)*3,(xnum-1)*(ynum-1)*3);    
    S = umfpack(L2,'\',R);
else
    L2 = sparse(I,J,L,(xnum-1)*(ynum-1)*3,(xnum-1)*(ynum-1)*3);
    S = L2\R;
end

% Compute the solution residual
res = R-L2*S;
resg = norm(res);

% Reload solution
vx=zeros(ynum+1,xnum);
vy=zeros(ynum,xnum+1);
pr=zeros(ynum-1,xnum-1);
% Initialize residual arrays
resx=zeros(ynum+1,xnum);Rx=resx;
resy=zeros(ynum,xnum+1);Ry=resy;
resc=zeros(ynum-1,xnum-1);Rc=resc;

for i=1:1:ynum-1
    for j=1:1:xnum-1
        % Indexes for P,vx,vy
        ivx=((j-1)*(ynum-1)+(i-1))*3+1;
        ivy=ivx+1;
        ipr=ivx+2;
        % Reload Vx
        if (j<xnum-1)
            vx(i+1,j+1)=S(ivx);
            resx(i+1,j+1)=res(ivx);
            Rx(i+1,j+1)=R(ivx);
        end
        % Reload Vy
        if (i<ynum-1)
            vy(i+1,j+1)=S(ivy);
            resy(i+1,j+1)=res(ivy);
            Ry(i,j)=R(ivy);
        end
        % Reload P
        pr(i,j)=S(ipr)*pscale;
        if pr(i,j)<0;
            pr(i,j) = prfirst(2);
        end
        % Continuity residual
        resc(i,j)=res(ipr);
        Rc(i,j)=R(ipr);
    end
end

resx = norm(resx);
resy = norm(resy);
resc = norm(resc); % Deviation from zero (i.e. mass conservation)

% Velocity boundary conditions

% Apply vx boundary conditions
% Left,Right Boundary
vx(:,1) = bleft(:,1) + bleft(:,2).*vx(:,2);
vx(:,xnum) = bright(:,1) + bright(:,2).*vx(:,xnum-1);
% Top, Bottom Boundary
vx(1,:) = btop(1:xnum,1)' + btop(1:xnum,2)'.*vx(2,:);
vx(ynum+1,:) = bbottom(1:xnum,1)' + bbottom(1:xnum,2)'.*vx(ynum,:);

% Apply vy boundary conditions ynum,xnum+1
% Left,Right Boundary
vy(:,1) = bleft(1:ynum,3) + bleft(1:ynum,4).*vy(:,2);
vy(:,xnum+1) = bright(1:ynum,3) + bright(1:ynum,4).*vy(:,xnum);
% Top, Bottom Boundary
vy(1,:) = btop(:,3)' + btop(:,4)'.*vy(2,:);
vy(ynum,:) = bbottom(:,3)' + bbottom(:,4)'.*vy(ynum-1,:);
