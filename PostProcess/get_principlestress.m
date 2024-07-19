function sii = get_principlestress(G)

for i = 2:G.ynum
    for j = 2:G.xnum
        
                sii(i-1,j-1)=sqrt(G.sxx1(i-1,j-1)^2+(G.sxy1(i-1,j-1)^2+G.sxy1(i,j-1)^2+G.sxy1(i-1,j)^2+G.sxy1(i,j)^2)/4);
                
    end
end