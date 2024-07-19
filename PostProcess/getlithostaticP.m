function [x,y,pr1] = getlithostaticP(G,bleft,bright,btop,bbottom,bintern)

% Make sure this is the same for your model if calculating overpressure
gx = 0;
gy=9.81;

prfirst(1)=0;
% prfirst(2) = boundary condition value
prfirst(2)=1e5;  % Pa  (for sticky layer condition set to 1e5)

% No need to use
useSuiteSparse = false;

yr2sec = 365.25*24*3600;


% Set very high viscosity, no movement, effectively a static model should
% yield vx and vy = 0.
etas0 = 1e30+zeros(size(G.etas1));
etan0 = 1e30+zeros(size(G.etan1));

gridx = G.gridx;
gridy = G.gridy;
xnum = G.xnum;
ynum = G.ynum;

% Computing right part of mechanical viscoelastic equations
% x-Stokes
RX1=zeros(ynum+1,xnum);
% y-Stokes
RY1=zeros(ynum,xnum+1);
% continuity
RC1=zeros(ynum-1,xnum-1);

% Grid points cycle
% No elastic stresses, purely viscous
for i=2:ynum
    for j=2:xnum

        % Right part of x-Stokes Equation
        if(j<xnum)
            RX1(i,j)=-gx*(G.rho1(i,j)+G.rho1(i-1,j))/2;
        end

        % Right part of y-Stokes Equation
        if(i<ynum)
            RY1(i,j)=-gy*(G.rho1(i,j)+G.rho1(i,j-1))/2;
        end
    end
end

% Solving of Stokes and Continuity equations on nodes
% and computing residuals
% by calling function Stokes_Continuity_solver_grid()
% with viscoelastic numerical viscosity
% and modified right parts
[vx1,resx1,vy1,resy1,pr1,resc1]=newFast_Stokes_Continuity_solver_sandbox2(prfirst,etas0,etan0,xnum,ynum,gridx,gridy,RX1,RY1,RC1,bleft,bright,btop,bbottom,bintern,useSuiteSparse);
% [G.vx1,resx1,G.vy1,resy1,G.pr1,resc1]=Fast_Stokes_Continuity_solver_sandbox(prfirst,etas0,etan0,xnum,ynum,gridx,gridy,RX1,RY1,RC1,bleft,bright,btop,bbottom,bintern);
disp(['Stokes Relative Residuals: x, y, cont: ',num2str([max(resx1(:)) max(resy1(:)) max(resc1(:))])])
disp(['Maximum velocity cm/yr: vx, vy: ',num2str(max(abs(vx1(:)))*100*yr2sec),' ',num2str(max(abs(vy1(:)))*100*yr2sec)])

% Caclulate center grid points for pressure nodes
x = (G.gridx(1:end-1) + G.gridx(2:end))/2;
y = (G.gridy(1:end-1) + G.gridy(2:end))/2;