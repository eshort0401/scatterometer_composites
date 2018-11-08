function disOut=plotTimeSeries(tranTimeSeries)

close all;

plotLegend=0;
savePlots=1;
sigLevel=0.5;

cMap=cbrewer('qual','Dark2',8,'spline');

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

figure('units','centimeters','pos',[0 0 6.5 5.5]);
hold on;
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + 2*ti(1);
bottom = outerpos(2) + 2*ti(2);
ax_width = outerpos(3) - 2.25*ti(1) - 2.25*ti(3);
ax_height = outerpos(4) - 2.375*ti(2) - 2.375*ti(4);
ax.Position = [left bottom ax_width ax_height];

distances(1)=1;
[~, disMinIndex]=min(abs(tranTimeSeries.distance-200));
distances(2)=disMinIndex;
[~, disMinIndex]=min(abs(tranTimeSeries.distance-400));
if tranTimeSeries.pProj(disMinIndex)<sigLevel
    distances=[distances disMinIndex];
end
% if max(tranTimeSeries.distance)>600
%     [~, disMinIndex]=min(abs(tranTimeSeries.distance-600));
%     if tranTimeSeries.pProj(disMinIndex)<sigLevel
%         distances=[distances disMinIndex];
%     end
% end

legendInfo=cell(1,length(distances));

% Find first statistically insignificant p-Value.
% k=find(tranTimeSeries.pProj>sigLevel,1,'first');
% 
% if ~isempty(k) || k~=1
%     distances(3)=k-2;
% else
%     distances=length(tranTimeSeries.distance);
% end

% Plot dotted line to indicate 0
plot([0 48],[0 0],'--','Color',[0 0 0],'LineWidth',1);

for i=1:length(distances)
    
    interpX = 0:0.1:48;  
    interpY = interp1([tranTimeSeries.times-24 tranTimeSeries.times ...
        tranTimeSeries.times+24 tranTimeSeries.times+48],...
        [tranTimeSeries.timeSeries(distances(i),:)...
        tranTimeSeries.timeSeries(distances(i),:) ... 
        tranTimeSeries.timeSeries(distances(i),:) ...
        tranTimeSeries.timeSeries(distances(i),:)],...
        interpX,'spline');
    
%     h(i)=plot(interpX,interpY,'Color',cMap(i,:),'LineWidth',1,'LineStyle','-');
    h(i)=plot(0:0.1:48,tranTimeSeries.harmFun{distances(i)}(0:0.1:48),...
        'Color',cMap(i,:),'LineWidth',1,'LineStyle','-');
    plot([tranTimeSeries.times tranTimeSeries.times+24],...
        [tranTimeSeries.timeSeries(distances(i),:) tranTimeSeries.timeSeries(distances(i),:)],...
        'LineStyle','none','Marker','o','MarkerEdgeColor',cMap(i,:),'MarkerSize',4);
    
    legendInfo{i}=sprintf('%6.2f km',tranTimeSeries.distance(distances(i)));
    
end

axis([0 48 -2 2]);
yticks(-2:.5:2);
xlim([0,48]);
xticks(0:4:48);

set(gca,'FontSize',12,'FontName','Times New Roman');
xlabel('LST','FontSize',12,'FontName','Times New Roman');
ylabel('m/s','FontSize',12,'FontName','Times New Roman');

if plotLegend
    legend(h(1:end),legendInfo,'FontSize',10,'FontName','Times New Roman');
end

if savePlots
    print(gcf,'-dsvg',[folderLoc,'/plot_',tranTimeSeries.label]);
end

disOut=tranTimeSeries.distance(distances);
    
end

