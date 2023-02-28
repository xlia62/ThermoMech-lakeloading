function fid = out2fastscape(input_dir,outfile,sframe,eframe,xlims)

yr2sec = 365.25*24*3600;

nframes = length(sframe:eframe);

c = 0;
for iframe = sframe:eframe
    c = c + 1;
    % Load the marker file for the iframe
    load([input_dir,'/grids_',num2str(iframe),'.mat'],'gridt','timesum')
    t(c) = timesum/yr2sec/1e6;  % time stamp in Myrs
    ix = gridt(1,:)/1000;     % convert to km
    iux = gridt(4,:)*yr2sec;  % convert to m/yr
    iuy = gridt(5,:)*yr2sec;  % convert to m/yr
    
    % interpolate according to imposed limits, keeping x position constant for
    % all subsquent time steps.
    if c == 1
        k = ix>=xlims(1) & ix<=xlims(2);
        x = ix(k);
        nx = length(x);
        ux = zeros(nframes,nx);
        uy = ux;
        ux(c,:) = iux(k);
        uy(c,:) = iuy(k);
    else
        ux(c,:) = interp1(ix,iux,x);
        uy(c,:) = interp1(ix,iuy,x);
    end
    
end

x = x - x(1);
uy = -uy;
if t(end) == t(end-1)
    x(end) = [];
    t(end) = [];
    ux(end) = [];
    uy(end) = [];
end
fid = writeFastScapeVelocity([input_dir,'/',outfile],x,t,ux,uy);



