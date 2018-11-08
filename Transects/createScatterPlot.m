function [perts,meanWinds,b] = createScatterPlot(tranCell)

close all;

perts=nan(length(tranCell),1);
meanWinds=nan(length(tranCell),1);
icons=['+','o','*','x','s','d','^','v','>','<','p','h']; 

hold on;

cMap=lines(3);

for i=1:8
    
    perts(i)=max(abs(tranCell{i}{1}.pertProj));
    % Average both times of day across transect to get mean wind
    meanWinds(i)=mean([tranCell{i}{1}.compProj, tranCell{i}{2}.compProj]);
    plot(meanWinds(i),perts(i),icons(mod(i-1,8)+1),'LineWidth',1.25,...
        'MarkerEdgeColor',cMap(1,:),'MarkerFaceColor',cMap(1,:));
    
end

for i=9:16
    
    perts(i)=max(abs(tranCell{i}{1}.pertProj));
    % Average both times of day across transect to get mean wind
    meanWinds(i)=mean([tranCell{i}{1}.compProj, tranCell{i}{2}.compProj]);
    plot(meanWinds(i),perts(i),icons(mod(i-1,8)+1),'LineWidth',1.25,...
        'MarkerEdgeColor',cMap(2,:),'MarkerFaceColor',cMap(2,:));
    
end

for i=17:24
    
    perts(i)=max(abs(tranCell{i}{1}.pertProj));
    % Average both times of day across transect to get mean wind
    meanWinds(i)=mean([tranCell{i}{1}.compProj, tranCell{i}{2}.compProj]);
    plot(meanWinds(i),perts(i),icons(mod(i-1,8)+1),'LineWidth',1.25,...
        'MarkerEdgeColor',cMap(3,:),'MarkerFaceColor',cMap(3,:));

end

xlim([-8 8]);
ylim([0 2.5]);
yticks(0:.5:2.5);

mW=[ones(length(meanWinds),1),meanWinds];
b=mW\perts;
x=min(meanWinds):0.01:max(meanWinds);

plot(x,b(1)+b(2)*x,'--k');

end

