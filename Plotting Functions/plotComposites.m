function handles=plotComposites(bin,plotArgs,startLat,endLat,startLon,endLon,folderLoc)

close all;

if nargin<=2
    % Use full grid
    startLon=bin.x(1);
    endLon=bin.x(end);
    startLat=bin.y(1);
    endLat=bin.y(end);
end

if nargin==1 || ~isstruct(plotArgs)
    clear plotArgs;
    plotArgs.composite=0;
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
    plotArgs.figureX=15; % In cm
    plotArgs.figureY=15;
   
    plotArgs.savePlots=0;
    plotArgs.labels=1;
    plotArgs.format='-dpng';
   
    plotArgs.compThresh=0.2;
    plotArgs.Rthresh=.2;
    % Skip this many grid cells before plotting an arrow.
    plotArgs.quiverSubSample=8;
    plotArgs.alphaFDR=0.05;
    plotArgs.simpleSig=0;
    plotArgs.tick=10;
    
end
    
% Set folder to save plots
    
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

%--------------------------------------------------------------------------
% Produce the plots specified by the input options.
%--------------------------------------------------------------------------

% Determine step sizes.
dLon=bin.x(2)-bin.x(1);
dLat=bin.y(2)-bin.y(1);

if nargin>2
    startLon=startLon+dLon/2;
    startLat=startLat+dLat/2;
    endLon=endLon-dLon/2;
    endLat=endLat-dLat/2;
end

% Define subbins
if (startLon~= bin.x(1) || endLon~=bin.x(end) || ...
        startLat~=bin.y(1) || endLat~=bin.y(end))
    
    % Check subgrid is actually within original grid
    if (startLon<bin.x(1) || endLon>bin.x(end) || ...
        startLat<bin.y(1) || endLat>bin.y(end))
        error('Specified grid outside range of data in bin.');
    end
    
    fprintf('Determining subgrids. \n');
    
    ratio=(endLat-startLat)./(bin.x(end)-bin.x(1));
    plotArgs.quiverSubSample=ceil(plotArgs.quiverSubSample*ratio);
    plotArgs.tick=ceil(plotArgs.tick*ratio);
    
    clear ratio;

    bin=subBin(bin,startLat,endLat,startLon,endLon);
    
end

% Define meshgrid for plotting winds and heatMaps
[X, Y]=meshgrid(bin.x,bin.y);

if bin.binArgs.LST==1
    timeScale='LST';
else
    timeScale='UTC';
end

MJOstring='';

% Set start and end times for the purposes of plot labels and file names.
startTime=bin.dateCell{1};
endTime=bin.dateCell{end};

if plotArgs.savePlots
   currentTime=clock;
   folderName=[num2str(floor(startLon)) '_' num2str(startTime(1)) '_' num2str(startTime(2)) '_run_' ...
       num2str(currentTime(1)) '_' num2str(currentTime(2)) '_' num2str(currentTime(3)) '_'...
       num2str(currentTime(4)) '_' num2str(currentTime(5)) '_' num2str(floor(currentTime(6))) '_' bin.binArgs.dataSource];
   mkdir(folderLoc,folderName);
end

m_proj('equidistant cylindrical','longitudes',[startLon endLon], ...
'latitudes',[startLat endLat]);

% Determine where there is insufficient data for coherent plot.
maxCount=max(bin.countComp(:));
insuffData=bin.countComp./maxCount<plotArgs.compThresh;

% Only plot where there is sufficient data.
[~ , ~, posCountTimes]=ind2sub(size(bin.countComp),find(~insuffData));
posCountTimes=unique(posCountTimes);

% Determine statistical significance
if plotArgs.simpleSig
    statSig=bin.pValue<=plotArgs.alphaFDR;
    statSig=repmat(statSig,[1 1 bin.NtPerDay]);
else
    [statSig, pStar]=FDRsignificance(bin.pValue,plotArgs.alphaFDR);
    statSig=repmat(statSig,[1 1 bin.NtPerDay]);
end
    
close all;

handlesN=1;

fprintf('Plotting. \n');

% Produce composite plots if required.
if plotArgs.composite
    
    bin.uComp(insuffData | bin.uComp==0 | isinf(bin.uComp))=NaN;
    bin.vComp(insuffData | bin.uComp==0 | isinf(bin.uComp))=NaN;

    compSpeed=sqrt(bin.uComp.^2+bin.vComp.^2);
    axisMin=floor(min(compSpeed(:)));
    axisMax=ceil(max(compSpeed(:)));
    
%     if axisMax-axisMin<6
%         step=.5;
%     else
%         step=1;
%     end
    
    axisMax=8;
    axisMin=0;
    step=1;

    for j=1:length(posCountTimes)
        
        i=posCountTimes(j);
        handles(handlesN)=plotWinds(X,Y,bin.uComp(:,:,i),bin.vComp(:,:,i),axisMin,step,axisMax,compSpeed(:,:,i),plotArgs);
        handlesN=handlesN+1;
        [plotTitle, figFileName]=makeTitles(i,startTime,endTime,bin.binArgs.dataSource,timeScale,'composite',bin.dt);
        title(plotTitle,'FontSize',12,'FontWeight','Normal','FontName','Times New Roman');
        if plotArgs.savePlots
%             saveas(gcf,[folderLoc,folderName,'/',figFileName],'fig');
            print(gcf,plotArgs.format,[folderLoc,folderName,'/',figFileName],'-r200');
            close(gcf);
        end
    end
    clear compSpeed
end

if plotArgs.pComp
    
    bin.pComp(isinf(bin.pComp))=NaN;
    
    axisMin=.25;
    axisMax=2;
    step=0.5;

    for j=1:length(posCountTimes)
        i=posCountTimes(j);
        handles(handlesN)=contourPlotPrecip(X,Y,bin.pComp(:,:,i),axisMin,axisMax,step,plotArgs.tick,true,5,2);
        handlesN=handlesN+1;
        [plotTitle, figFileName]=makeTitles(i,startTime,endTime,bin.binArgs.dataSource,timeScale,'rain rate',bin.dt);
        title(plotTitle,'FontSize',12,'FontWeight','Normal','FontName','Times New Roman');
        if plotArgs.savePlots
            print(gcf,plotArgs.format,[folderLoc,folderName,'/',figFileName],'-r200');
            close(gcf);
        end
    end
    
    clear compSpeed

end

if plotArgs.pPertComp
    
    bin.pPertComp(isinf(bin.pPertComp))=NaN;
    
    axisMin=-1;
    axisMax=1;
    step=0.1;

    for j=1:length(posCountTimes)
        i=posCountTimes(j);
        handles(handlesN)=contourPlot(X,Y,bin.pPertComp(:,:,i),axisMin,axisMax,step,plotArgs.tick,true,4);
        handlesN=handlesN+1;
        [plotTitle, figFileName]=makeTitles(i,startTime,endTime,bin.binArgs.dataSource,timeScale,'rain rate',bin.dt);
        title(plotTitle,'FontSize',12,'FontWeight','Normal','FontName','Times New Roman');
        if plotArgs.savePlots
            print(gcf,plotArgs.format,[folderLoc,folderName,'/',figFileName],'-r100');
            close(gcf);
        end
    end
    
    clear compSpeed

end

% Produce composite plots if required.
if plotArgs.mean
    
    bin.uComp(insuffData | bin.uComp==0 | isinf(bin.uComp))=NaN;
    bin.vComp(insuffData | bin.uComp==0 | isinf(bin.uComp))=NaN;
    
    uMean=nanmean(bin.uComp,3);
    vMean=nanmean(bin.vComp,3);

    compSpeed=sqrt(uMean.^2+vMean.^2);
    axisMin=floor(min(compSpeed(:)));
    axisMax=ceil(max(compSpeed(:)));
    
    axisMin=0;
    axisMax=8;
    
    if axisMax-axisMin<6
        step=.5;
    else
        step=1;
    end
   
    handles(handlesN)=plotWinds(X,Y,uMean,vMean,axisMin,step,axisMax,compSpeed,plotArgs);
    handlesN=handlesN+1;
    [plotTitle, figFileName]=makeTitles(1,startTime,endTime,bin.binArgs.dataSource,timeScale,'daily_composite',bin.dt);
    if plotArgs.labels
        title(plotTitle,'FontSize',12,'FontWeight','Normal','FontName','Times New Roman');
    end
    if plotArgs.savePlots
        print(gcf,plotArgs.format,[folderLoc,folderName,'/',figFileName],'-r200');
        close(gcf);
    end
    
    clear compSpeed

end

if plotArgs.pertComp
        
    bin.uPertComp(insuffData | bin.uPertComp==0 | isinf(bin.uPertComp) | ~statSig)=NaN;
    bin.vPertComp(insuffData | bin.vPertComp==0 | isinf(bin.vPertComp) | ~statSig)=NaN;
    
    pertCompSpeed=sqrt(bin.uPertComp.^2+bin.vPertComp.^2);
    axisMin=floor(min(pertCompSpeed(:)));
    axisMax=ceil(max(pertCompSpeed(:)));
    
    if axisMax-axisMin<8
        step=.5;
    elseif axisMax-axisMin<16
        step=1;
    else
        step=2;
    end
    
    axisMax=4;
    axisMin=0;
    step=.5;
    
    for j=1:length(posCountTimes)
        
        i=posCountTimes(j);
        handles(handlesN)=plotWinds(X,Y,bin.uPertComp(:,:,i),bin.vPertComp(:,:,i),axisMin,step,axisMax,pertCompSpeed(:,:,i),plotArgs);
        handlesN=handlesN+1;
        [plotTitle, figFileName]=makeTitles(i,startTime,endTime,bin.binArgs.dataSource,timeScale,'perturbation_composite',bin.dt);
        if plotArgs.labels
            title(plotTitle,'FontSize',12,'FontWeight','Normal','FontName','Times New Roman');
        end
        if plotArgs.savePlots
            print(gcf,plotArgs.format,[folderLoc,folderName,'/',figFileName],'-r200');
            close(gcf);
        end
        
    end
    
    clear pertCompSpeed;

end

% Produce heat maps if required.
if plotArgs.countMap
   
    axisMin=floor(min(min(bin.countComp(bin.countComp>0)))/20)*20;
    axisMax=ceil(max(max(bin.countComp(bin.countComp>0)))/20)*20;
    step=20;
% 
%     axisMin=0;
%     axisMax=3;
%     step=1;

   
    for j=1:length(posCountTimes)

        i=posCountTimes(j);
        handles(handlesN)=contourPlot(X,Y,bin.countComp(:,:,i),axisMin,axisMax,step,plotArgs.tick,1);
        handlesN=handlesN+1;
        [plotTitle, figFileName]=makeTitles(i,startTime,endTime,bin.binArgs.dataSource,timeScale,'composite count',bin.dt);
%         title(plotTitle);
        if plotArgs.savePlots
            saveas(gcf,[folderLoc,folderName,'/',figFileName],'fig');
            print(gcf,plotArgs.format,[folderLoc,folderName,'/',figFileName],'-r200');
            close(gcf);
        end

    end
end

% Produce stdev plot if required.
if plotArgs.uVar || plotArgs.vVar
    bin.uPertVar(insuffData | bin.uPertVar<=0 | isinf(bin.uPertVar))=NaN;
    bin.uPertVar=sqrt(bin.uPertVar);

    bin.uPertVar(insuffData | bin.vPertVar<=0 | isinf(bin.vPertVar))=NaN;
    bin.vPertVar=sqrt(bin.vPertVar);

    axisMin=floor(min([bin.vPertVar(:); bin.uPertVar(:)]));
    axisMax=ceil(max([bin.vPertVar(:); bin.uPertVar(:)]));
end

if plotArgs.uVar
    
%     if axisMax-axisMin<8
%         step=.5;
%     elseif axisMax-axisMin<16
%         step=1;
%     elseif axisMax-axisMin<64
%         step=4;
%     else
%         step=8;
%     end

    axisMin=0;
    step=.5;
    axisMax=4;
    
    for j=1:length(posCountTimes)

        i=posCountTimes(j);
        handles(handlesN)=contourPlot(X,Y,bin.uPertVar(:,:,i),axisMin,axisMax,step,plotArgs.tick,1);
        handlesN=handlesN+1;
        [plotTitle, figFileName]=makeTitles(i,startTime,endTime,bin.binArgs.dataSource,timeScale,'u standard deviation',bin.dt);
%         title(plotTitle);
        if plotArgs.savePlots
            saveas(gcf,[folderLoc,folderName,'/',figFileName],'fig');
            print(gcf,plotArgs.format,[folderLoc,folderName,'/',figFileName],'-r200');
            close(gcf);
        end

    end

end

% Produce variance plot if required.
if plotArgs.vVar
        
%     if axisMax-axisMin<8
%         step=.5;
%     elseif axisMax-axisMin<16
%         step=1;
%     elseif axisMax-axisMin<64
%         step=4;
%     else
%         step=8;
%     end
%     

    axisMin=0;
    step=.5;
    axisMax=4;

    for j=1:length(posCountTimes)

        i=posCountTimes(j);
        contourPlot(X,Y,bin.vPertVar(:,:,i),axisMin,axisMax,step,plotArgs.tick,1);
        [plotTitle, figFileName]=makeTitles(i,startTime,endTime,bin.binArgs.dataSource,timeScale,'v standard deviation',bin.dt);
%         title(plotTitle);
        if plotArgs.savePlots
            saveas(gcf,[folderLoc,folderName,'/',figFileName],'fig');
            print(gcf,plotArgs.format,[folderLoc,folderName,'/',figFileName],'-r800');
            close(gcf);
        end

    end

end

if plotArgs.coVar

    bin.pertCoVar(insuffData | bin.pertCoVar==0 | isinf(bin.pertCoVar))=NaN;
    bin.pertCoVar=bin.pertCoVar./(bin.uPertVar.*bin.vPertVar);
    
    axisMin=-1;
    axisMax=1;
    
    step=.25;
        
    for j=1:length(posCountTimes)

        i=posCountTimes(j);
        contourPlot(X,Y,bin.pertCoVar(:,:,i),axisMin,axisMax,step,plotArgs.tick,1,1);
        [plotTitle, figFileName]=makeTitles(i,startTime,endTime,bin.binArgs.dataSource,timeScale,'correlation',bin.dt);
%         title(plotTitle);
        if plotArgs.savePlots
            saveas(gcf,[folderLoc,folderName,'/',figFileName],'fig');
            print(gcf,plotArgs.format,[folderLoc,folderName,'/',figFileName],'-r800');
            close(gcf);
        end

    end

end

if plotArgs.tVar

    bin.tVar(insuffData | bin.tVar==0 | isinf(bin.tVar))=NaN;
    bin.tVar=sqrt(bin.tVar);
    
    axisMin=0;
    axisMax=ceil(max(bin.tVar(:)*60)/5)*5;
    
    step=5;
        
    for j=1:length(posCountTimes)

        i=posCountTimes(j);
        contourPlot(X,Y,bin.tVar(:,:,i)*60,axisMin,axisMax,step,plotArgs.tick,'time std deviation');
        [plotTitle, figFileName]=makeTitles(i,startTime,endTime,bin.binArgs.dataSource,timeScale,'standard deviation',bin.dt);
%         title(plotTitle);
        if plotArgs.savePlots
            saveas(gcf,[folderLoc,folderName,'/',figFileName],'fig');
            print(gcf,plotArgs.format,[folderLoc,folderName,'/',figFileName],'-r800');
            close(gcf);
        end

    end

end

if plotArgs.tComp
    
    bin.tComp(insuffData | bin.tComp==0 | isinf(bin.tComp))=NaN;
    axisMin=0;
    axisMax=ceil(bin.dt./(60));
    step=10;
    
    for j=1:length(posCountTimes)
        
        i=posCountTimes(j);
        minutesSinceBinStart=mod(bin.tComp(:,:,i)-(startTime(4)+startTime(5)/60+startTime(6)/3600 + (i-1)*bin.dt./(3600)),24)*60;
        contourPlot(X,Y,minutesSinceBinStart,axisMin,axisMax,step,plotArgs.tick,1,1,0);
        [plotTitle, figFileName]=makeTitles(i,startTime,endTime,bin.binArgs.dataSource,timeScale,'average time',bin.dt);
%         title(plotTitle);
        if plotArgs.savePlots
            saveas(gcf,[folderLoc,folderName,'/',figFileName],'fig');
            print(gcf,plotArgs.format,[folderLoc,folderName,'/',figFileName],'-r800');
            close(gcf);
        end

    end

end

if plotArgs.div
    
    % Note divergence values are much less than compositing accuracy.
    % Satellite tracks appear visible in divergence fields.
    bin.div(insuffData | isinf(bin.div))=NaN;
    step=10^(-4);
    axisMin=floor(min(bin.div(:))/step)*step;
    axisMax=ceil(max(bin.div(:))/step)*step;
    
    for j=1:length(posCountTimes)
        
        i=posCountTimes(j);
        contourPlot(X,Y,bin.div(:,:,i),axisMin,axisMax,step,plotArgs.tick,' divergence ',1);
        [plotTitle, figFileName]=makeTitles(i,startTime,endTime,bin.binArgs.dataSource,timeScale,' divergence ',bin.dt);
%         title(plotTitle);
        if plotArgs.savePlots
            saveas(gcf,[folderLoc,folderName,'/',figFileName],'fig');
            print(gcf,plotArgs.format,[folderLoc,folderName,'/',figFileName],'-r800');
            close(gcf);
        end

    end

end

if plotArgs.ellipseMean
    
    time=9;
    
    uMean=bin.diurnal.uCoef(:,:,2).*cos((2*pi/24)*time)+bin.diurnal.uCoef(:,:,3).*sin((2*pi/24)*time);
    vMean=bin.diurnal.vCoef(:,:,2).*cos((2*pi/24)*time)+bin.diurnal.vCoef(:,:,3).*sin((2*pi/24)*time);
    compSpeed=sqrt(uMean.^2+vMean.^2);
    
    goodFit=(bin.diurnal.uR>=plotArgs.Rthresh & bin.diurnal.vR>=plotArgs.Rthresh);
    
    uMean(~goodFit | all(insuffData,3))=NaN;
    vMean(~goodFit | all(insuffData,3))=NaN;
    compSpeed(~goodFit | all(insuffData,3))=NaN;
    
    axisMin=floor(min(compSpeed(:)));
    axisMax=ceil(max(compSpeed(:)));
    
    if axisMax-axisMin<=2
        step=.25;
    elseif axisMax-axisMin<6
        step=.5;
    else
        step=1;
    end

    plotWinds(X,Y,uMean,vMean,axisMin,step,axisMax,compSpeed,plotArgs.quiverSubSample,plotArgs.tick);
    [plotTitle, figFileName]=makeTitles(1,startTime,endTime,bin.binArgs.dataSource,timeScale,'Ellipse mean',bin.dt);
%     title(plotTitle,'FontSize',12,'FontWeight','Normal','FontName','Times New Roman');
    if plotArgs.savePlots
        saveas(gcf,[folderLoc,folderName,'/',figFileName],'fig');
        print(gcf,plotArgs.format,[folderLoc,folderName,'/',figFileName],'-r800');
        close(gcf);
    end

    clear compSpeed

end

if plotArgs.tMax
    
    axisMin=0;
    axisMax=12;
    step=1;        
    goodFit=(bin.diurnal.uR>=plotArgs.Rthresh & bin.diurnal.vR>=plotArgs.Rthresh);
    bin.diurnal.tMax(~goodFit | all(insuffData,3))=NaN;
    
    contourPlot(X,Y,bin.diurnal.tMax,axisMin,axisMax,step,plotArgs.tick,' time max perturbation ',2,1);
    [plotTitle, figFileName]=makeTitles(1,startTime,endTime,bin.binArgs.dataSource,timeScale,' time max perturbation ',bin.dt);
    title(plotTitle);
    if plotArgs.savePlots
        saveas(gcf,[folderLoc,folderName,'/',figFileName],'fig');
        print(gcf,plotArgs.format,[folderLoc,folderName,'/',figFileName],'-r800');
        close(gcf);
    end

end

end
