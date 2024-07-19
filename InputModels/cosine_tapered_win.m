function y=cosine_tapered_win(n,r)

m = n*r/2;

j = 0:n-1;

y = ones(size(j));

k = j <= m-1;
y(k) = 0.5*(1-cos(pi*j(k)/m));

k = j >= n-m;
y(k) = 0.5*(1-cos(pi*(n-j(k)-1)/m));

end