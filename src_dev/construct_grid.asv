function G = construct_grid(G)

% This function constructs the grid structure according to the grid
% structure

G.cellnum = (G.xnum-1)*(G.ynum-1);  % Number of cells

switch G.grid_type
            
    case 'boundary'  % Increase grid density near boundaries
        
        % Set boundary node spacing
        dx = 6e3;  % meters
        dy = 6e3;  % meters
        G.xb = 10;   % 10 nodes for xboundary
        G.yb = 10;   % 10 nodes for yboundary
                
        x1 = zeros(1,G.xb);
        y1 = zeros(1,G.yb);
        for i = 2:G.xb
            x1(i) = x1(i-1) + dx;
        end
        
        for i = 2:G.yb
            y1(i) = y1(i-1) + dy;
        end
        
        x3 = G.xsize - x1(end:-1:1);
        y3 = G.ysize - y1(end:-1:1);       
        
        x2 = linspace(x1(end),x3(1),G.xnum-2*G.xb + 2);
        y2 = linspace(y1(end),y3(1),G.ynum-2*G.yb + 2);
                
        G.x = [x1 x2(2:end-1) x3];
        G.y = [y1 y2(2:end-1) y3];
       
    case 'nonuniform'
        
        [G.x,G.y,G.xstp,G.ystp] = build_grid2(G.xnum,G.ynum,G.xsize,G.ysize,G.bx,G.by,G.wx,G.wy,G.f);
        G.x = G.x'; G.y = G.y';
        
    case 'nonuniform_balance'
        
        [G.x,G.y,G.xstp,G.ystp] = build_grid3(G.xnum,G.ynum,G.xsize,G.ysize,G.bx,G.by,G.wx,G.wy,G.f);
        G.x = G.x'; G.y = G.y';
        
    case 'nonuniform_open'
        
        [G.x,G.y,G.xstp,G.ystp] = build_grid2(G.xnum,G.ynum,G.xsize,G.ysize,G.bx,G.by,G.wx,G.wy,G.f);        
        G.x = G.x'; G.y = G.y';
        
    case 'specified'
        
    case 'adaptive'
        
    case 'gerya_subduction'
        
        [G.x,G.y,G.xstp,G.ystp] = build_grid4(G.xnum,G.ynum,G.xsize,G.ysize,G.bx,G.by,G.wx,G.wy,G.f,G.xpos);        
        G.x = G.x'; G.y = G.y';
        
        
    otherwise % default uniform
        
        G.grid_type = 'uniform';
        
        % Create vectors for nodal points positions (basic nodes)
        G.x = linspace(0,G.xsize,G.xnum ); % Horizontal
        G.y = linspace(0,G.ysize,G.ynum ); % Vertical        
        
end

% Grid spacing
G.dx = diff(G.x);
G.dy = diff(G.y);

% Node connection defining the pressure centered cell
[G.X,G.Y]= meshgrid(G.x,G.y);
X = G.X(:);
Y = G.Y(:);
G.conxy = connect_quad(G.ynum,G.xnum);
G.xx = X(G.conxy);
G.yy = Y(G.conxy);

% Cell centers
G.xc = G.x(1:end-1) + G.dx/2;
G.yc = G.y(1:end-1) + G.dy/2;


