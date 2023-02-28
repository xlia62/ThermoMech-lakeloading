% Boundary conditions file used to apply velocity and thermal boundary 
% conditions to a model. This script can be edited to apply sepecialized
% boundary conditions.

% To be updated to include stress boundary option

% Staggered Grid
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

% -------------------------------------------------------------------------
%% Boundary Conditions -- Parameters
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Pressure boundary conditions
% -------------------------------------------------------------------------
% prfirst(1) = boundary condition mode:
% 0 - pressure in one cell definition
% 1 - pressure at the top and in the bottom of the channel
prfirst(1)=0;
% prfirst(2) = boundary condition value
prfirst(2)=1e5;  % Pa  (for sticky layer condition set to 1e5)

% Moving grid requires special mass conserving boundary condition
moving_grid = true;
regridx = false;

% If moving grid, update the bottom temperature condition according to
% adiabat?
update_tbottom = false;

% Fixed Internal boundary condtion
internalbound = false;

% Moving Internal boundary condition
movebound = false;

% Use external boundary condition as a permeable lower boundary condition
external_velobound = false;
external_tempbound = false;
zext = 1000e3;   % Where to apply external condition, must be >> gridy(end)
   
% -------------------------------------------------------------------------
% Remove horizontal boundary conditions after period of time remove_time%
% -------------------------------------------------------------------------
remove_time = false;%3.5;   % Ma, set to false to turn-off% Strechtching try 3.5Ma%false

%------------------------------------------------------------------------------
% Horizontal Velocity Boundary (xleft, xright, xbottom, xtop) (cm/yr)
% Can be a function or a scalar -- Edit this section for your conditions

% Set the vx velocity to v0 be uplied to the upper d0 km, and then use a
% linear ramp down to 0 cm/yr from d0 km to d0 + d1 km.
global v0
if isempty(v0)
    v0 = 0.15;%1.5;
end

vxleft = -v0; 
vxright = v0;

vxtop = 0;
vxbot = 0;

% Vertical Velocity Boundary (yleft, yright, ybottom, top) (cm/yr), 
% NOTE: negative velocity means upward flow towards surface
vyleft = 0; 
vyright = 0;
vytop = 0;
vybot = 0;

%------------------------------------------------------------------------------
% Type of boundary condition
% Fixed (Rigid) use true; for freeslip = false
xtopfixed = false;
ytopfixed = true;
xbotfixed = false;
ybotfixed = true;

xleftfixed = true;
yleftfixed = true;
xrightfixed = true;
yrightfixed = true;

% Special conditions
% Upper, Lower boundaries: Free slip + Prescribed inward velocity (vertical shortening)
% Velocities are determined according to equations by Laoi and Gerya,
% Tectonophysics, 2014

% External temperature condition, in this implimentation moving_grid must
% be false
% Lower Boundary external fixed (grids can't change)
if external_tempbound && moving_grid
   error('External Temperature Boundary Condition is not compatible with a moving grid')
end

% Moving grid is a special condition and needs to ensure conservation of
% mass
if moving_grid
    ytopfixed = true;
    ybotfixed = true;
    vytop = sticky_layer*(vxright-vxleft)/G.xsize;
    vybot = (vxright-vxleft)*(sticky_layer-G.ysize)/G.xsize;
end

% External boundary is akin to permeable (though not set through a stress
% boundary condition)
if external_velobound
   xbotfixed = false;
   dz = gridy(ynum)-gridy(ynum-1);
   vybot = 0;
end

% -------------------------------------------------------------------------
%% Boundary Conditions -- Velocity
% -------------------------------------------------------------------------

% Convert to SI units
vxleft = vxleft/100/yr2sec;
vyleft = vyleft/100/yr2sec;
vxright = vxright/100/yr2sec;
vyright = vyright/100/yr2sec;

vxtop = vxtop/100/yr2sec;
vytop = vytop/100/yr2sec;
vxbot = vxbot/100/yr2sec;
vybot = vybot/100/yr2sec;

%---------------------------------------------
% Initialize
%---------------------------------------------
btop = zeros(xnum+1,4);
bbottom = zeros(xnum+1,4);
bleft = zeros(ynum+1,4);
bright = zeros(ynum+1,4);

%---------------------------------------------
% Upper boundary
%---------------------------------------------
% vx(1,j)=btop(j,1)+vx(2,j)*btop(j,2)
% vy(1,j)=btop(j,3)+vy(2,j)*btop(j,4)

% Default is free-slip
btop(:,2) = 1;
btop(:,4) = 1;

if xtopfixed
    btop(:,1) = vxtop;
    btop(:,2) = 0;
end

if ytopfixed
    btop(:,3) = vytop;
    btop(:,4) = 0;
end   

%---------------------------------------------
% Lower boundary
%---------------------------------------------
% vx(ynum+1,j)=bbottom(j,1)+vx(ynum,j)*bbottom(j,2)
% vy(ynum,j)=bbottom(j,3)+vy(ynum-1,j)*bbottom(j,4)

% Default is free-slip
bbottom(:,2)= 1;
bbottom(:,4)= 1;

if xbotfixed
    bbottom(:,1) = vxbot;
    bbottom(:,2) = 0;
end

if ybotfixed
    bbottom(:,3) = vybot;
    bbottom(:,4) = 0;
end

if external_velobound
    bbottom(:,3) = 0;
    bbotoom(:,4) = zext/(zext+dz);
end

%---------------------------------------------
% Left boundary: 
%---------------------------------------------
% vx(i,1)=bleft(i,1)+vx(i,2)*bleft(i,2)
% vy(i,1)=bleft(i,3)+vy(i,2)*bleft(i,4)

% Default is free-slip
bleft(:,2) = 1;
bleft(:,4) = 1;

if xleftfixed
    bleft(:,1) = vxleft;
    bleft(:,2) = 0;
end

if yleftfixed
    bleft(:,3) = vyleft;
    bleft(:,4) = 0;
end

%---------------------------------------------
% Right boundary
%---------------------------------------------
% vx(i,xnum)=bright(i,1)+vx(i,xnum-1)*bbright(i,2)
% vy(i,xnum+1)=bright(i,3)+vx(i,xnum)*bbright(i,4)

% Default is free-slip
bright(:,2) = 1;
bright(:,4) = 1;

if xrightfixed
    bright(:,1) = vxright;
    bright(:,2) = 0;
end

if yrightfixed
    bright(:,3) = vyright;
    bright(:,4) = 0;
end

% -------------------------------------------------------------------------
%% Boundary Conditions -- Thermal
% -------------------------------------------------------------------------

%---------------------------------------------
% Initialize
%---------------------------------------------
btopt = zeros(xnum,2);
bbottomt = btopt;
bleftt = zeros(ynum,2);
brightt = bleftt;

%---------------------------------------------
% Upper Boundary fixed
%---------------------------------------------
% tk(1,j)=btopt(j,1)+tk(2,j)*btop(j,2)
btopt(:,1) = ttop;
btopt(:,2) = 0;

%---------------------------------------------
% Lower Boundary fixed
%---------------------------------------------
% tk(ynum,j)=bbottomt(j,1)+tk(ynum-1,j)*bbottomt(j,2)
bbottomt(:,1) = tbottom;
bbottomt(:,2) = 0;

if external_tempbound
    dz = gridy(ynum)-gridy(ynum-1);
    Text = tpotential + tgrad*(zext+gridy(ynum)-sticky_layer);
    bbottomt(:,1) = Text*dz/(zext+dz);
    bbottomt(:,2) = zext/(zext+dz);
end

%---------------------------------------------
% Left, Right boundaries: symmetry (zero flux)
%---------------------------------------------
for i=1:1:ynum
    % Left boundary
    % tk(i,1)=bleftt(i,1)+bleftt(i,2)*tk(i,2);
    bleftt(i,1)=0;
    bleftt(i,2)=1;
    % Right boundary
    % tk(i,xnum)=brightt(i,1)+brightt(i,2)*tk(i,xnum-1);
    brightt(i,1)=0;
    brightt(i,2)=1;    
end

% -------------------------------------------------------------------------
%% Boundary Condition -- Internal
% -------------------------------------------------------------------------

% Internal boundary condition: prescribed velocity of "mobile wall"
% Postion are described in terms of nodes
bintern = zeros(8,4);
bintern(1,:)=-1;      % Horizontal position of vx nodes with prescrbed velocity (no susch condition if negative)
bintern(2,:)=0;       % Min vertical position
bintern(3,:)=0;       % Max vertical position
bintern(4,:)=-0;      % Prescribed shortening velocity, m/s
bintern(5,:)=-1;      % Horizontal position of vy nodes with prescrbed velocity (no susch condition if negative)
bintern(6,:)=0;       % Min vertical position
bintern(7,:)=0;       % Max vertical position
bintern(8,:)=0;       % Prescribed vertical velocity, m/s
