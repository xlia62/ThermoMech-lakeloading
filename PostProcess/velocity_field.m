% Visualize the flow field
function velocity_field(filename, timestep,x,y)


% Load grid file
% run to timestep
tmp = [filename, '/grids_',num2str(timestep),'.mat'];
load(tmp)

animate_matprop(filename,timestep,timestep,500,500,markercolors(16),[],[min(x) max(x)],[min(y) max(y)],false)

% meshgrid (km)
[X,Y] = meshgrid(x,y);
[vx,vy] = velocityplot(X*1000,Y*1000,G);

hold on
quiver(X,Y,vx,vy,2,'k')
hold off


yr2sec = 365.25*24*3600;
Vx = vx*100*yr2sec;
Vy = vy*100*yr2sec;

%figure
%pcolor(X,Y,sqrt(Vx.^2 + Vy.^2))
%shading interp
%colorbar
