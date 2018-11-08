function totalData = plotErrorHist(edges,time,varargin)
%ERRORHISTOGRAMS Summary of this function goes here
%   Detailed explanation goes here

% Set figure window
figure('units','centimeters','pos',[0 0 14 10]);

cMap=cbrewer('qual','Set1',9,'pchip');

counts=nan(length(edges)-1,length(varargin));
counts(:,1)=histcounts(varargin{1}(:,:,1),edges)';

totalData=sum(sum(~isnan(varargin{1}(:,:,1)),1),2);

for i=2:length(varargin)
    counts(:,i)=histcounts(varargin{i}(:,:,time),edges)';
end

step=edges(2)-edges(1);
labels=edges(1)+step/2:step:edges(end)-step/2;
errorBar=bar(labels,counts);

for i=1:length(varargin)
    errorBar(i).FaceColor = cMap(i,:);
    errorBar(i).LineWidth = 0.25;
end

set(gca,'FontSize',12,'FontName','Times New Roman');
set(gca,'XTick',edges,'FontSize',12,'FontName','Times New Roman')


ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + 1.75*ti(1);
bottom = outerpos(2) + 2*ti(2);
ax_width = outerpos(3) - 2.25*ti(1) - 2.25*ti(3);
ax_height = outerpos(4) - 2.375*ti(2) - 2.375*ti(4);
ax.Position = [left bottom ax_width ax_height];
