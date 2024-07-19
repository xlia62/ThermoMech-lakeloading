function [mx,my,mxn,myn,mcxn,mcyn,mi,meta,mtk] = reflectmarkers(gridx,gridy,gridcx,gridcy,MX,MY,MI,META,MTK)

% This function moves markers from the right side of model space to the
% left side

mx = [];
my = [];
mxn = [];
myn = [];
mcxn = [];
mcyn = [];
mi = [];
meta = [];
mtk = [];
k = MX > gridx(end);
newmarkers = sum(k);

if newmarkers
    
    % New x-position of the markers on the left side
    mx = MX(k) - gridx(end);
    
    % Copy rest of information to new markers
    my = MY(k);
    mi = MI(k);
    meta = META(k);
    mtk = MTK(k);
    
    mxn = zeros(size(mx));
    myn = mxn;
    mcxn = mxn;
    mcyn = mxn;
    
    xnum = length(gridx);
    ynum = length(gridy);
    
    for mm1 = 1:newmarkers
        
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
            xn=double(int16((xnmax+xnmin)./2-0.5));
            if(gridx(xn)>mx(mm1))
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
        % Save horizontal index
        mxn(mm1)=xn;
        
        % Now define the index for the central nodes (Presure, SXX, etc..)
        if (mx(mm1) < gridcx(xn+1))
            xn = xn-1;
        end
        if (xn<1)
            xn=1;
        end
        if (xn>xnum-2)
            xn=xnum-2;
        end
        mcxn(mm1) = xn;
        
        % Find vertical index
        ynmin=1;
        ynmax=ynum;
        while ((ynmax-ynmin)>1)
            % !!! SUBTRACT 0.5 since int16(0.5)=1
            yn=double(int16((ynmax+ynmin)./2-0.5));
            if(gridy(yn)>my(mm1))
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
        % Save Vertical index
        myn(mm1)=yn;
        
        % Now define the index for the central nodes (Presure, SXX, etc..)
        if (my(mm1) < gridcy(yn+1))
            yn = yn-1;
        end
        if(yn<1)
            yn = 1;
        end
        if(yn>ynum-2)
            yn = ynum-2;
        end
        mcyn(mm1) = yn;
        
    end
    
end

