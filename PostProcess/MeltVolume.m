
% Model A melt volume.
c = 1;
[tA(c),VmeltA(c)]=totalmeltextracted('../MCRplumeA',20,500,500,[420 580],[0 400],3000,2900);

for i = [40:20:760 769]
    c = c + 1;
    [tA(c),VmeltA(c)]=totalmeltextracted('../MCRplumeA',i,500,500,[420 580],[0 400],3000,2900);
    if VmeltA(c) < VmeltA(c-1)
        VmeltA(c) = VmeltA(c-1);
    end
end

% Model B melt volume.
c = 1;
[tB(c),VmeltB(c)]=totalmeltextracted('../MCRplumeB',20,500,500,[420 580],[0 400],3000,2900);

for i = [40:20:920 934]
    c = c + 1;
    [tB(c),VmeltB(c)]=totalmeltextracted('../MCRplumeB',i,500,500,[420 580],[0 400],3000,2900);
    if VmeltB(c) < VmeltB(c-1)
        VmeltB(c) = VmeltB(c-1);
    end
end

% Model C melt volume.
c = 1;
[tC(c),VmeltC(c)]=totalmeltextracted('../MCRplumeC',20,500,500,[420 580],[0 400],3000,2900);

for i = [40:20:1000]
    c = c + 1;
    [tC(c),VmeltC(c)]=totalmeltextracted('../MCRplumeC',i,500,500,[420 580],[0 400],3000,2900);
    if VmeltC(c) < VmeltC(c-1)
        VmeltC(c) = VmeltC(c-1);
    end
end

% Model D melt volume.
c = 1;
[tD(c),VmeltD(c)]=totalmeltextracted('../MCRplumeD',20,500,500,[420 580],[0 400],3000,2900);

for i = [40:20:320 330]
    c = c + 1;
    [tD(c),VmeltD(c)]=totalmeltextracted('../MCRplumeD',i,500,500,[420 580],[0 400],3000,2900);
    if VmeltD(c) < VmeltD(c-1)
        VmeltD(c) = VmeltD(c-1);
    end
end

% Model E melt volume.
c = 1;
[tE(c),VmeltE(c)]=totalmeltextracted('../MCRplumeE',20,500,500,[420 580],[0 400],3000,2900);

for i = [40:20:540 549]
    c = c + 1;
    [tE(c),VmeltE(c)]=totalmeltextracted('../MCRplumeE',i,500,500,[420 580],[0 400],3000,2900);
    if VmeltE(c) < VmeltE(c-1)
        VmeltE(c) = VmeltE(c-1);
    end
end

% Model F melt volume.
c = 1;
[tF(c),VmeltF(c)]=totalmeltextracted('../MCRplumeF',20,500,500,[420 580],[0 400],3000,2900);

for i = [40:20:560 567]
    c = c + 1;
    [tF(c),VmeltF(c)]=totalmeltextracted('../MCRplumeF',i,500,500,[420 580],[0 400],3000,2900);
    if VmeltF(c) < VmeltF(c-1)
        VmeltF(c) = VmeltF(c-1);
    end
end

figure
plot(tA,VmeltA,'linewidth',2)
hold on
plot(tB,VmeltB,'linewidth',2)
plot(tC,VmeltC,'linewidth',2)
plot(tD,VmeltD,'linewidth',2)
plot(tE,VmeltE,'linewidth',2)
plot(tF,VmeltF,'linewidth',2)

xlabel('Model Time (My)')
ylabel('Melt Volume (km^{2}/km)')
title('Modeled Extracted Cumulative Melt Volume')
set(gca,'FontSize',14)

xlim([0 10])
ylim([0 1325])

legend('A','B','C','D','E','F')
grid

nt = 100;
t = linspace(0,10,nt);
Vmelt = zeros(6,nt);
t(end) = 10.1;

Vmelt(1,:) = interp1(tA,VmeltA,t);
Vmelt(2,:) = interp1(tB,VmeltB,t);
Vmelt(3,:) = interp1(tC,VmeltC,t);
Vmelt(4,:) = interp1(tD,VmeltD,t);
Vmelt(5,:) = interp1(tE,VmeltE,t);
Vmelt(6,:) = interp1(tF,VmeltF,t);

[C,h] = contour(repmat(t,6,1),Vmelt,repmat([1:6]',1,nt),[1.01 2 3 4.01 5 5.99],'linewidth',2);
th = clabel(C,h,'manual','FontSize',14);
set(th,'VerticalAlignment','bottom')

st = 'ABCDEF';
cmp = ([1 0 1;1 0.84 0;1 0 0;0 0.5 0;0 0 1;0 0 0]);
for i = 1:length(th)
    k = round(str2double(th(i).String));
    th(i).String = st(k);
    th(i).Color = cmp(k,:);
end


xlim([0 10])
ylim([0 1325])
caxis([1 6])
colormap(cmp)
hold on
plot([0 10],[1050 1050],'k--','linewidth',2) % Estimated Volume of flood basalt
xlabel('Model Time (My)')
ylabel('Volume (km^{3}/km)')
title('Modeled Extracted Cumulative Melt Volume')
set(gca,'FontSize',14)
grid


