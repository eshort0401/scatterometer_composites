function [pStarActive, pStarInactive, pStarT1, pStarT2]=plotMJOanalysis(binAvActive,binAvInactive)

savePlots=1;
alpha=0.05;

plotArgs.composite=1;
plotArgs.mean=0;
plotArgs.pertComp=1;
plotArgs.countMap=0;
plotArgs.pComp=0;
plotArgs.pPertComp=0;
plotArgs.uVar=0;
plotArgs.vVar=0;
plotArgs.coVar=0;
plotArgs.tComp=0;
plotArgs.tVar=0;
plotArgs.div=0;
plotArgs.ellipseMean=0;
plotArgs.tMax=0;
plotArgs.figureX=13; % In cm
plotArgs.figureY=6;

plotArgs.savePlots=0;
plotArgs.labels=0;
plotArgs.format='-dpng';

plotArgs.compThresh=0.2;
plotArgs.Rthresh=.4;
% Skip this many grid cells before plotting an arrow.
plotArgs.quiverSubSample=8;
plotArgs.alphaFDR=0.05;
plotArgs.simpleSig=0;
plotArgs.tick=10;

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

if savePlots
   currentTime=clock;
   folderName=[num2str(currentTime(1)) '_' num2str(currentTime(2)) '_' num2str(currentTime(3)) '_'...
       num2str(currentTime(4)) '_' num2str(currentTime(5)) '_' num2str(floor(currentTime(6)))];
   mkdir(folderLoc,folderName);
end

% Calculate statistically significant difference between active and
% inactive for both T1 and T2.
pValueT1=hotelling(binAvActive,binAvInactive,1,1);
pValueT2=hotelling(binAvActive,binAvInactive,9,9);

pValueActive=binAvActive.pValue;
pValueInactive=binAvInactive.pValue;

[~, pStarActive]=FDRsignificance(pValueActive,alpha);
[~, pStarInactive]=FDRsignificance(pValueInactive,alpha);

pValueT1(pValueActive>pStarActive | pValueInactive>pStarInactive)=nan;
pValueT2(pValueActive>pStarActive | pValueInactive>pStarInactive)=nan;

[~, pStarT1]=FDRsignificance(pValueT1,alpha);
[~, pStarT2]=FDRsignificance(pValueT2,alpha);

[i1, j1]=find(pValueT1<=pStarT1);
[i2, j2]=find(pValueT2<=pStarT2);

% Produce active phase plots
hActive=plotComposites(binAvActive,plotArgs);
set(0, 'currentfigure', hActive(3));
m_plot(binAvActive.x(i1),binAvActive.y(j1),'o','markers',2,'Color',[.1,.1,0.9]);

set(0, 'currentfigure', hActive(4));
m_plot(binAvActive.x(i2),binAvActive.y(j2),'o','markers',2,'Color',[.1,.1,0.9]);

if savePlots
    print(hActive(1),'-dpng',[folderLoc,folderName,'/a1comp.png'],'-r200');
    print(hActive(2),'-dpng',[folderLoc,folderName,'/a2comp.png'],'-r200');
    print(hActive(3),'-dpng',[folderLoc,folderName,'/a1pert.png'],'-r200');
    print(hActive(4),'-dpng',[folderLoc,folderName,'/a2pert.png'],'-r200');
end

% Define meshgrid for plotting winds and heatMaps
[X, Y]=meshgrid(binAvActive.x,bin.y);

% Produce inactive phase plots
hInactive=plotComposites(binAvInactive,plotArgs);
set(0, 'currentfigure', hInactive(3));
m_plot(binAvInactive.x(i1),binAvInactive.y(j1),'o','markers',2,'Color',[.1,.1,0.9]);
% contourPlot(binAvActive.x,binAvActive.y,z,axisMin,axisMax,step,tick,labels,div,pcolour)

set(0, 'currentfigure', hInactive(4));
% m_plot(binAvInactive.x(i2),binAvInactive.y(j2),'o','markers',2,'Color',[.1,.1,0.9]);

if savePlots
    print(hInactive(1),'-dpng',[folderLoc,folderName,'/i1comp.png'],'-r200');
    print(hInactive(2),'-dpng',[folderLoc,folderName,'/i2comp.png'],'-r200');
    print(hInactive(3),'-dpng',[folderLoc,folderName,'/i1pert.png'],'-r200');
    print(hInactive(4),'-dpng',[folderLoc,folderName,'/i2pert.png'],'-r200');
end

m_plot(binAvInactive.x(i1),binAvInactive.y(j1),'o','markers',2,'Color',[.1,.1,0.9]);


% Produce mean plot

% Plot difference between MJO active and inactive
plotArgs.pertComp=0;
plotArgs.composite=0;
plotArgs.mean=1;

% binAvDiff=binAvActive;
% binAvDiff.uComp=binAvActive.uComp-binAvInactive.uComp;
% binAvDiff.vComp=binAvActive.vComp-binAvInactive.vComp;

% hDiff=plotComposites(binAvDiff,plotArgs);

% if savePlots
%     print(hDiff(1),'-dpng',[folderLoc,folderName,'/diffMean.png'],'-r200');
% end

end

