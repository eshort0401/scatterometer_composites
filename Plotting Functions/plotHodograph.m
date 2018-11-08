function plotHodograph(u,v,c,t)

close all

figure('units','centimeters','pos',[0 0 15 8]);
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + 1.5*ti(1);
bottom = outerpos(2) + 1.5*ti(2);
ax_width = outerpos(3) - 2.5*ti(1) - 2.5*ti(3);
ax_height = outerpos(4) - 2*ti(2) - 2*ti(4);
ax.Position = [left bottom ax_width ax_height];

uMinusC=u-c(1);
vMinusC=v-c(2);

plot(squeeze(uMinusC(33,25,:,t)), squeeze(vMinusC(33,25,:,t)),'-or');
hold on;
plot(c(1),c(2),'-xr');
plot(squeeze(u(33,25,:,t)),squeeze(v(33,25,:,t)),'-ok');
axis([-25 25 -25 0])
plot([0 0],[0 -25],'--','color',[.4 .4 .4]);
% plot([25 -25],[0 0],'--','color',[.4 .4 .4]);
set(gca,'FontSize',12','FontName','Times New Roman')

end

