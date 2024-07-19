function [mx,my,mxn,myn,mcxn,mcyn,nmark] = fillemptycell(gridx,gridy,gridcx,gridcy,MXN,MYN,mxcell,mycell,filllastrow)

xnum = length(gridx);
ynum = length(gridy);

markspercell = mxcell*mycell;
maxmarks = 20*markspercell;
mx = zeros(maxmarks,1);
my = zeros(maxmarks,1);
mxn = zeros(maxmarks,1);
myn = zeros(maxmarks,1);
nmark = 1;

if filllastrow  % Fill only the last row
    
    for j = 1:xnum-1
        
        % Find all markers associated with the x-grid line
        k1 = MXN == j;
        myk1 = MYN(k1);
        
        i = ynum-1;
        
        k2 = myk1 == i;
        
        if isempty(k2) || sum(k2) == 0
            dx = gridx(j+1) - gridx(j);
            dy = gridy(i+1) - gridy(i);
            
            [x,y] = quassidistribute(dx,dy,mxcell,mycell,2);
            
            if nmark > maxmarks
                mx = [mx;mx];
                my = [my;my];
                mxn = [mxn;mxn];
                myn = [myn;myn];
                maxmarks = maxmarks + maxmarks;
            end
            
            mx(nmark:nmark+markspercell-1) = gridx(j) + x;
            my(nmark:nmark+markspercell-1) = gridy(i) + y;
            mxn(nmark:nmark+markspercell-1) = j;
            myn(nmark:nmark+markspercell-1) = i;
            
            nmark = nmark + markspercell;
            
        end
        
    end
    
    
else
    
    for j = 1:xnum-1
        
        % Find all markers associated with the x-grid line
        k1 = MXN == j;
        myk1 = MYN(k1);
        
        for i = 1:ynum-1
            
            k2 = myk1 == i;
            
            if isempty(k2) || sum(k2) == 0
                dx = gridx(j+1) - gridx(j);
                dy = gridy(i+1) - gridy(i);
                
                [x,y] = quassidistribute(dx,dy,mxcell,mycell,2);
                
                if nmark > maxmarks
                    mx = [mx;mx];
                    my = [my;my];
                    mxn = [mxn;mxn];
                    myn = [myn;myn];
                    maxmarks = maxmarks + maxmarks;
                end
                
                mx(nmark:nmark+markspercell-1) = gridx(j) + x;
                my(nmark:nmark+markspercell-1) = gridy(i) + y;
                mxn(nmark:nmark+markspercell-1) = j;
                myn(nmark:nmark+markspercell-1) = i;
                
                nmark = nmark + markspercell;
                
            end
            
        end
    end
    
end

nmark = nmark - 1;
mx = mx(1:nmark);
my = my(1:nmark);
mxn = mxn(1:nmark);
myn = myn(1:nmark);
mcxn = zeros(size(mxn));
mcyn = zeros(size(myn));

for mm1 = 1:nmark
    
    xn = mxn(mm1);
    yn = myn(mm1);
    
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


