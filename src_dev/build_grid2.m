function [gridx,gridy,xstp,ystp] = build_grid2(xnum,ynum,xsize,ysize,bx,by,wx,wy,f)

% Assamble a grid with the following properties:
%
% 1) A horizontal grid that has a centered fine grid spanning width wx of
% model space and is bounded by a non-uniform grid of equal size on each
% side:
%           
% x0        x1               x2        xsize
% +---------+ <---- wx ----> +---------+
% |         |                |         |
%

%
% 2) A vertical grid
%

% x0 = x(1) grid position (m)
% b = desired uniform grid spacing in fine part of the grid (m)
% wx = width of fine uniform grid (m)

% Defining average initial gridsteps
xstp=xsize./(xnum-1);
ystp=ysize./(ynum-1);

% Defining gridline positions for irregular basic grid

%% Horizontal grid
gridx=zeros(xnum,1);

% Number of xnodes required for fine grid
nx = fix(wx/bx);

if nx > xnum
    error('Number of x-nodes required for fine grid is equal to x-grid size; increase the number of x-nodes')
end
% Center node
xc = fix((xnum+1)*f);

% Initial position of the area
x1 = xc - (nx-mod(nx,2))/2;
if x1 < 1
    error('Increase horizontal fine grid spacing, or increase number of nodes, or decrease the width')
end
x2 = x1 + nx;

gridx(x1) = xsize*f - wx/2;

% Assemble fine grid
for i = x1+1 : x2
    gridx(i) = gridx(i-1) + bx;
end

% Define factor of grid spacing increase to the right
% of high resolution area
D = xsize - gridx(x2); % distance to be covered by non-uniform grid

N = xnum-x2; % number of grid steps to be used in the grid

if N == 0 || wx == 0
    % We have a uniform grid
    gridx = linspace(0,xsize,xnum)'; 
else
    % Iterative search of F
    F=1.1;
    for i=1:1:100
        F=(1+D./bx.*(1-1/F)).^(1/N);
    end
    % Define position of nodal points
    for i=x2:xnum
        gridx(i)=gridx(i-1)+bx*F.^(i-x2);
    end
    gridx(xnum)=xsize;
        
    % Define factor of grid spacing increase to the left
    % of high resolution area
    
    D=gridx(x1); % distance to be covered by non-uniform grid
    N=x1-1; % number of grid steps to be used in the grid
    
    % Iterative search of F
    F=1.1;
    for i=1:100
        F=(1+D./bx.*(1-1/F)).^(1/N);
    end
    
    % Define position of nodal points
    gridx(1)=0;
    for i=2:N
        gridx(i)=gridx(i-1)+bx*F.^(N+2-i);
    end
    
end
%% Vertical grid
gridy=zeros(ynum,1);

% Number of xnodes required for fine grid
y1 = fix(wy/by);

if y1 > ynum
    error('Number of y-nodes required for fine grid is equal to y-grid size; increase the number of x\y-nodes')
end

% Define regular step in high resolution area
for i=2:y1
    gridy(i)=gridy(i-1)+by;
end

% Define factor of grid spacing increase from the bottom
% of high resolution area
D=ysize-gridy(y1); % distance to be covered by non-uniform grid
N = ynum - y1 ; % number of grid steps to be used in the grid

if N == 0 || wy == 0
    % We have a uniform grid
    gridy = linspace(0,ysize,ynum )'; 
else
    % Iterative search of F
    F=1.1;
    for i=1:1:100
        F=(1+D./by.*(1-1/F)).^(1/N);
    end
    
    % Define position of nodal points
    for i=y1+1:ynum
        gridy(i)=gridy(i-1)+by*F.^(i-y1);
    end
    %    gridy(ynum)=ysize;
end