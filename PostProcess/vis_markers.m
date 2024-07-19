function [px,py,markprop] = vis_markers(gridx,gridy,npx,npy,MX,MY,MPROP)
%UNTITLED2 Summary of this function goes here
% npx,npy = number of image pixels in x and y direction
%   Detailed explanation goes here

ngrid = 2;
% Visualizing marker type

% Pixel grid resolution
px = linspace(gridx(1),gridx(end),npx);
py = linspace(gridy(1),gridy(end),npy);
sxstp = px(2) - px(1);
systp = py(2) - py(1);

% Process markers
markprop = NaN+ones(npy,npx);
markdis = 1e20*ones(npy,npx);
% loop over markerrs
marknum = length(MX);

for mm1 = 1:1:marknum
    
    % Define pixel cell
    m1=fix(MX(mm1)/sxstp)+1;
    m2=fix(MY(mm1)/systp)+1;
    if (m1<1)
        m1=1;
    end
    if (m1>npx-1)
        m1=npx-1;
    end
    if (m2<1)
        m2=1;
    end
    if (m2>npy-1)
        m2=npy-1;
    end
    % Define indexes of surrounding pixels
    m10min=m1-ngrid;
    if (m10min<1)
        m10min=1;
    end
    m10max=m1+1+ngrid;
    if (m10max>npx)
        m10max=npx;
    end
    m20min=m2-ngrid;
    if (m20min<1)
        m20min=1;
    end
    m20max=m2+1+ngrid;
    if (m20max>npy)
        m20max=npy;
    end
    
    mx = MX(mm1);
    my = MY(mm1);
    % Update pixels around the marker
    for m10 = m10min:1:m10max
        for m20 = m20min:1:m20max
            % Check distance to current pixel
            dx=mx-(m10-1)*sxstp;
            dy=my-(m20-1)*systp;
            dd=sqrt(dx*dx+dy*dy);
            if(dd<markdis(m20,m10))
                markprop(m20,m10)=MPROP(mm1);
                markdis(m20,m10)=dd;
            end
        end
    end
end

end

