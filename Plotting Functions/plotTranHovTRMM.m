function plotTranHovTRMM(tranCell,times)

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
time=[times-24 times times+24];

for i=1:length(tranCell)
    tranAll(i,:)=tranCell{i}.compProj;
end

tranAll=vertcat(tranAll,tranAll,tranAll);

% axisMin=-2;
% axisStep=.2;
% axisMax=2;

scale=.5.^(4:-1:0);

[X, Y]=meshgrid(distance,time);

cMap=cbrewer('seq','Greys',length(scale)+1,'pchip');
cMap=cMap(2:end,:);

figureScale=max(tranCell{1}.distance)/995.45;
% figure('units','centimeters');
figure('units','centimeters','pos',[0 0 8 8]);
ax = gca;

%     outerpos = ax.OuterPosition;
%     ti = ax.TightInset; 
%     left = outerpos(1) + 1.5*ti(1);
%     bottom = outerpos(2) + 1.5*ti(2);
%     ax_width = outerpos(3) - 2.5*ti(1) - 2.5*ti(3);
%     ax_height = outerpos(4) - 2*ti(2) - 2*ti(4);
%     ax.Position = [left bottom ax_width ax_height];

[C, h]=contour(X,Y,log2(tranAll),log2(scale),'LineWidth',1.5);
colormap(cMap);
c=colorbar;
caxis([log2(scale(1))-0.5, log2(scale(end))+0.5]);
set(c,'YTick',log2(scale),'YTickLabel',scale,'FontSize',12,'FontName','Times New Roman')
ylabel(c,'mm/h','FontSize',12,'FontName','Times New Roman');
% clabel(C,h,'FontSize',10,'FontName','Times New Roman','LabelSpacing',5000);
axis([0 max(tranCell{1}.distance) 0 48]);
ax.PlotBoxAspectRatio=[figureScale 1 1];
% set(ch,'edgecolor','none');
clear ch;

k=find(tranCell{i}.pProj>sigLevel,1,'first');
hold on;
if ~isempty(k) && k>1
    plot([tranCell{1}.distance(k) tranCell{1}.distance(k)],[0 48],...
        '--','Color',[0 0 0],'LineWidth',1);
end

% colormap(cMap);
% c=colorbar;
% caxis([axisMin,axisMax]);
% set(c,'YTick',axisMin:axisStep:axisMax,'FontSize',12,'FontName','Times New Roman')
% ylabel(c,'m/s','FontSize',12,'FontName','Times New Roman');
yticks(0:4:48);
xticks(0:200:1000);

set(gca,'FontSize',12','FontName','Times New Roman')

if savePlots
    print(gcf,'-dsvg',[folderLoc,'/hov_',tranCell{1}.label,'.svg']);
end
    
end

