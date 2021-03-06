function plotTranHov(tranCell,times)

% close all;

savePlots=1;
sigLevel=0.05;

% Personal Macbook.
if ismac
    folderLoc='/Volumes/Ewan''s Hard Drive/Figures/';
end

% Uni Unix box machines.
if isunix && not(ismac)
    username=char(java.lang.System.getProperty('user.name'));
    folderLoc=['/media/' username '/Ewan''s Hard Drive/Figures/'];
    clear username;
end

tranAll=zeros(length(tranCell),length(tranCell{1}.distance));
distance=tranCell{1}.distance;
time=[times-24 times times+24 times+48];

for i=1:length(tranCell)
    tranAll(i,:)=tranCell{i}.pertProj;
end

tranAll=vertcat(tranAll,tranAll,tranAll, tranAll);

axisMin=-2.5;
axisStep=.25;
axisMax=2.5;

[X, Y]=meshgrid(distance,time);

cMap=flipud(cbrewer('div','RdBu',length(axisMin:axisStep:axisMax)-1,'pchip'));

figureScale=max(tranCell{1}.distance)/995.45;

figure('units','centimeters','pos',[0 0 8 8]);
ax = gca;

%     outerpos = ax.OuterPosition;
%     ti = ax.TightInset; 
%     left = outerpos(1) + 1.5*ti(1);
%     bottom = outerpos(2) + 1.5*ti(2);
%     ax_width = outerpos(3) - 2.5*ti(1) - 2.5*ti(3);
%     ax_height = outerpos(4) - 2*ti(2) - 2*ti(4);
%     ax.Position = [left bottom ax_width ax_height];

[~, ch]=contourf(X,Y,tranAll,axisMin:axisStep:axisMax); 
axis([0 max(tranCell{1}.distance) 0 48]);
ax.PlotBoxAspectRatio=[figureScale 1 1];
set(ch,'edgecolor','none');
clear ch;

k=find(tranCell{i}.pProj>sigLevel,1,'first');
hold on;
if ~isempty(k) && k>1
    plot([tranCell{1}.distance(k) tranCell{1}.distance(k)],[0 48],...
        '--','Color',[0 0 0],'LineWidth',1);
end

colormap(cMap);
c=colorbar;
caxis([axisMin,axisMax]);
set(c,'YTick',axisMin:axisStep:axisMax,'FontSize',12,'FontName','Times New Roman')
ylabel(c,'m/s','FontSize',12,'FontName','Times New Roman');
yticks(0:4:48);
xticks(0:200:1000);

set(gca,'FontSize',12','FontName','Times New Roman')

if savePlots
    print(gcf,'-dpng',[folderLoc,'/hov_',tranCell{1}.label],'-r200');
end
    
end

