function [compMagError, pertCompMagError, statSig] = compareData(bin1,bin2,startLat,endLat,startLon,endLon,plotArgs)

if nargin<=2
    % Use full grid
    startLon=bin1.x(1);
    endLon=bin1.x(end);
    startLat=bin1.y(1);
    endLat=bin1.y(end);
end

if nargin<=3
    % Use full grid
    startLon=bin1.x(1);
    endLon=bin1.x(end);
    startLat=bin1.y(1);
    endLat=bin1.y(end);
end

if nargin<=6
    % Use default plot arguments
    plotArgs.compError=1;
    plotArgs.meanError=0;
    plotArgs.pertError=1;
    plotArgs.savePlots=0;
    plotArgs.format='-dpng';
    plotArgs.compThresh1=0;
    plotArgs.compThresh2=0;
     
    plotArgs.quiverSubSample=64;
    plotArgs.alphaFDR=0.05;
    tick=10;    
end

% Determine step sizes.
dLon=bin1.x(2)-bin1.x(1);
dLat=bin1.y(2)-bin1.y(1);

% Check two data sets are consistent
if (any(bin1.x ~= bin2.x) || any(bin1.y ~= bin2.y))
    error('Inconsistent grids')
elseif (bin1.NtPerDay ~= bin2.NtPerDay)
    error('Inconsistent times');
end

% Define subbins
if (startLon~= bin1.x(1) || endLon~=bin1.x(end) || ...
        startLat~=bin1.y(1) || endLat~=bin1.y(end))
    
    % Check subgrid is actually within original grid
    if (startLon<bin1.x(1) || endLon>bin1.x(end) || ...
        startLat<bin1.y(1) || endLat>bin1.y(end))
        error('Specified grid outside range of data in bin.');
    end
    
    fprintf('Determining subgrids. \n');
    
    ratio=(endLat-startLat)./(bin1.x(end)-bin1.x(1));
    plotArgs.quiverSubSample=ceil(plotArgs.quiverSubSample*ratio);
    tick=ceil(tick*ratio);
    
    clear ratio;

    bin1=subBin(bin1,startLat,endLat,startLon,endLon);
    bin2=subBin(bin2,startLat,endLat,startLon,endLon);
    
end

% Calculate error in comp

% Remove values without common obs.
bin1.uComp(bin1.countComp==0 | bin2.countComp==0)=0;
bin1.vComp(bin1.countComp==0 | bin2.countComp==0)=0;
bin2.uComp(bin1.countComp==0 | bin2.countComp==0)=0;
bin2.vComp(bin1.countComp==0 | bin2.countComp==0)=0;

compMagOne=sqrt(bin1.uComp.^2+bin1.vComp.^2);
compMagTwo=sqrt(bin2.uComp.^2+bin2.vComp.^2);

compMagError=compMagOne-compMagTwo;

% Calculate error in pertComp

% Remove values without common obs or statistical significance.
bin1.uPertComp(bin1.countComp==0 | bin2.countComp==0)=0;
bin1.vPertComp(bin1.countComp==0 | bin2.countComp==0)=0;
bin2.uPertComp(bin1.countComp==0 | bin2.countComp==0)=0;
bin2.vPertComp(bin1.countComp==0 | bin2.countComp==0)=0;

pertCompMagOne=sqrt(bin1.uPertComp.^2+bin1.vPertComp.^2);
pertCompMagTwo=sqrt(bin2.uPertComp.^2+bin2.vPertComp.^2);

pertCompMagError=pertCompMagOne-pertCompMagTwo;

%--------------------------------------------------------------------------
% Plot
%--------------------------------------------------------------------------
fprintf('Plotting. \n');

startTime=bin1.dateCell{1};

if plotArgs.savePlots
   currentTime=clock;
   folderName=['ERROR' num2str(startLon) '_' num2str(startTime(1)) '_' num2str(startTime(2)) '_run_' ...
       num2str(currentTime(1)) '_' num2str(currentTime(2)) '_' num2str(currentTime(3)) '_'...
       num2str(currentTime(4)) '_' num2str(currentTime(5)) '_' num2str(floor(currentTime(6)))];
   
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
   
    mkdir(folderLoc,folderName);
end

if bin1.binArgs.LST==1
    timeScale='LST';
else
    timeScale='UTC';
end

MJOstring=num2str(bin1.binArgs.mjoPhases(1));
for i=2:length(bin1.binArgs.mjoPhases)
    MJOstring=[MJOstring ', ' num2str(bin1.binArgs.mjoPhases(i))];
end

% Set start and end times for the purposes of plot labels and file names.
startTime=bin1.dateCell{1};
endTime=bin1.dateCell{end};

[X, Y]=meshgrid(bin1.x,bin1.y);
m_proj('equidistant cylindrical','longitudes',[bin1.x(1) bin1.x(end)], ...
'latitudes',[bin1.y(1) bin1.y(end)]);

% Only plot where there is a positive countComp.
[~ , ~, posCountTimes]=ind2sub(size(bin1.countComp),find(bin1.countComp>0 & bin2.countComp>0));
posCountTimes=unique(posCountTimes);

% Determine where there is insufficient data for coherent plot.
maxCount1=max(bin1.countComp(:));
maxCount2=max(bin2.countComp(:));

insuffData=(bin1.countComp./maxCount1<plotArgs.compThresh1 | ...
bin2.countComp./maxCount2<plotArgs.compThresh2);
    
close all;

% Produce composite error plots if required.
if plotArgs.compError
    
    compMagError(insuffData  | compMagError==0)=NaN;
    bin1.uComp(insuffData)=NaN;
    bin1.vComp(insuffData)=NaN;
    bin2.uComp(insuffData)=NaN;
    bin2.vComp(insuffData)=NaN;
    
%     if max(abs(compMagError(:)))<2
%         axisMin=-max(ceil(abs(compMagError(:)*5)))/5;
%         axisMax=max(ceil(abs(compMagError(:)*5)))/5;
%     else
%         axisMin=-max(ceil(abs(compMagError(:))));
%         axisMax=max(ceil(abs(compMagError(:))));
%     end
%         
%     if axisMax-axisMin<=4
%         step=.5;
%     else
%         step=1;
%     end

    axisMin=-0.5;
    axisMax=0.5;
    step=.1;
    
    for j=1:length(posCountTimes)
        i=posCountTimes(j);
        
        plotCompareWinds(X,Y,bin1.uComp(:,:,i),bin1.vComp(:,:,i),bin2.uComp(:,:,i),bin2.vComp(:,:,i),axisMin,step,axisMax,compMagError(:,:,i),plotArgs.quiverSubSample,tick);
        [plotTitle, figFileName]=makeTitlesComparison(i,startTime,endTime,bin1.binArgs.dataSource,bin2.binArgs.dataSource,timeScale,'difference',bin1.dt,MJOstring);
%         title(plotTitle,'FontSize',12,'FontWeight','Normal','FontName','Times New Roman');
        if plotArgs.savePlots
            saveas(gcf,[folderLoc,folderName,'/',figFileName],'fig');
            print(gcf,plotArgs.format,[folderLoc,folderName,'/',figFileName],'-r200');
            close(gcf);
        end

    end
    
end

% Produce mean wind error if required
if plotArgs.meanError
    
    % Calculate error in mean
    bin1.uComp(insuffData | bin1.uComp==0 | isinf(bin1.uComp))=NaN;
    bin1.vComp(insuffData | bin1.vComp==0 | isinf(bin1.uComp))=NaN;
    bin2.uComp(insuffData | bin2.uComp==0 | isinf(bin2.uComp))=NaN;
    bin2.vComp(insuffData | bin2.vComp==0 | isinf(bin2.uComp))=NaN;

    uMean1=nanmean(bin1.uComp,3);
    vMean1=nanmean(bin1.vComp,3);

    uMean2=nanmean(bin2.uComp,3);
    vMean2=nanmean(bin2.vComp,3);

    compSpeed1=sqrt(uMean1.^2+vMean1.^2);
    compSpeed2=sqrt(uMean2.^2+vMean2.^2);

    meanCompMagError=compSpeed1-compSpeed2;
    
%     if max(abs(compMagError(:)))<2
%         axisMin=-max(ceil(abs(compMagError(:)*5)))/5;
%         axisMax=max(ceil(abs(compMagError(:)*5)))/5;
%     else
%         axisMin=-max(ceil(abs(compMagError(:))));
%         axisMax=max(ceil(abs(compMagError(:))));
%     end
%         
%     if axisMax-axisMin<=4
%         step=.5;
%     else
%         step=1;
%     end

    axisMin=-10;
    axisMax=10;
    step=2;
        
    plotCompareWinds(X,Y,uMean1,vMean1,uMean2,vMean2,axisMin,step,axisMax,meanCompMagError,plotArgs.quiverSubSample,tick);
    [plotTitle, figFileName]=makeTitlesComparison(i,startTime,endTime,bin1.binArgs.dataSource,bin2.binArgs.dataSource,timeScale,'difference',bin1.dt,MJOstring);
%         title(plotTitle,'FontSize',12,'FontWeight','Normal','FontName','Times New Roman');
    if plotArgs.savePlots
        saveas(gcf,[folderLoc,folderName,'/',figFileName],'fig');
        print(gcf,plotArgs.format,[folderLoc,folderName,'/',figFileName],'-r200');
        close(gcf);
    end


    
end

% Produce perturbation composite error plots if required.
if plotArgs.pertError
    
    statSig1=FDRsignificance(bin1.pValue,plotArgs.alphaFDR);
    statSig1=repmat(statSig1,[1 1 bin2.NtPerDay]);
    statSig2=FDRsignificance(bin2.pValue,plotArgs.alphaFDR);
    statSig2=repmat(statSig2,[1 1 bin2.NtPerDay]);

    statSig=statSig1 & statSig2;
 
    clear statSig1 statSig2;
    
    pertCompMagError(insuffData  | ~statSig | pertCompMagError==0)=NaN;
    bin1.uPertComp(insuffData | ~statSig)=NaN;
    bin1.vPertComp(insuffData | ~statSig)=NaN;
    bin2.uPertComp(insuffData | ~statSig)=NaN;
    bin2.vPertComp(insuffData | ~statSig)=NaN;
    
%     if max(abs(pertCompMagError(:)))<=2
%         axisMin=-max(ceil(abs(pertCompMagError(:)*5)))/5;
%         axisMax=max(ceil(abs(pertCompMagError(:)*5)))/5;
%     else
%         axisMin=-max(ceil(abs(pertCompMagError(:))));
%         axisMax=max(ceil(abs(pertCompMagError(:))));
%     end
%     
%     if axisMax-axisMin<=4
%         step=.5;
%     else
%         step=1;
%     end

    axisMax=2;
    axisMin=-2;
    step=0.5;
    
    for j=1:length(posCountTimes)
        i=posCountTimes(j);
        
        plotCompareWinds(X,Y,bin1.uPertComp(:,:,i),bin1.vPertComp(:,:,i),bin2.uPertComp(:,:,i),bin2.vPertComp(:,:,i),axisMin,step,axisMax,pertCompMagError(:,:,i),plotArgs.quiverSubSample,tick);
        [plotTitle, figFileName]=makeTitlesComparison(i,startTime,endTime,bin1.binArgs.dataSource,bin2.binArgs.dataSource,timeScale,'difference in perturbations',bin1.dt,MJOstring);
%         title(plotTitle,'FontSize',12,'FontWeight','Normal','FontName','Times New Roman');
        if plotArgs.savePlots
            saveas(gcf,[folderLoc,folderName,'/',figFileName],'fig');
            print(gcf,plotArgs.format,[folderLoc,folderName,'/',figFileName],'-r200');
            close(gcf);   
        end

    end
    
end

end

