function [mxn,myn,cxn,cyn] = locatemarkeringrid(MX,MY,gridx,gridy,gridcx,gridcy,xnum,ynum)

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
    xn = floor((xnmax+xnmin)/2);
    %xn=double(int16((xnmax+xnmin)./2-0.5));
    if(gridx(xn)>MX)
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

mxn = xn;
% Now define the index for the central nodes (Presure, SXX, etc..)
if (MX < gridcx(xn+1))
    xn = xn-1;
end
if (xn<1)
    xn=1;
elseif (xn>xnum-2)
    xn=xnum-2;
end
    
cxn = xn;

ynmin=1;
ynmax=ynum;
while ((ynmax-ynmin)>1)
    % !!! SUBTRACT 0.5 since int16(0.5)=1
    % yn=double(int16((ynmax+ynmin)./2-0.5));
    yn = floor((ynmax+ynmin)/2);
    if(gridy(yn)>MY)
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

myn = yn;
% Now define the index for the central nodes (Presure, SXX, etc..)
if (MY < gridcy(yn+1))
    yn = yn-1;
end
if(yn<1)
    yn = 1;
elseif(yn>ynum-2)
    yn = ynum-2;
end

cyn = yn;
