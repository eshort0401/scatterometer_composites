function bin=binWRFprecip(fileNames,dLat,startLat,endLat,...
    dLon,startLon,endLon,NtPerDay,dateCell,ASCATcount,binArgs)
% plot_wind_dir_WRF.m ----------------------------------------------------------
% 
% Copyright Ewan Short 23/2/2016
%
% DESC: This function takes wind direction from netCDF format at time point in 
% column t and plots. WRF grid runs from -13.0142 to 6.0894 lat, 93.7142 to
% 152.8788 lon.
%
% INPUTS: 
% dLat and dLon are the
% size of lat and lon increments respectively, dt time increment in 
% seconds. Note the
% time variable in the netCDF data is measured in seconds since 1/1/1990
% midnight. tRange is a vector [startTime endTime] with components in
% same units. 
%
% OUTPUT: Output the binned u and v component winds as a complete data structure. 
%
%--------------------------------------------------------------------------

if nargin<11
    binArgs=struct;
    % Use LST (1) or UTC (0).
    binArgs.LST=1; 
    % Only bin certain MJO phases (1)
    binArgs.mjo=0;
    % Which phases?
    binArgs.mjoPhases=[1 2 3 4 5 6 7 8];
    binArgs.dataSource='WRF';
end

if nargin<10
    % Count corresponding to an ASCAT swath. Used for comparison purposes.
    ASCATcount=NaN;
end

close all;

% Determine number of time intervals. Check have been inputted correctly.
numTimes=length(dateCell);
numIntervals=numTimes/2;

if floor(numIntervals)~=numIntervals
    error('Odd number of date vectors. \n');
end

% Calculate total number of days and confirm an integer. 
bin.numDays=zeros(1,numIntervals);

for i=1:(numIntervals)
    bin.numDays(i)=etime(dateCell{2*i},dateCell{2*i-1})/(24*60*60);
    if floor(bin.numDays(i))~=bin.numDays
        error('Time intervals must be in whole days. \n');
    end
end

bin.numDaysTot=sum(bin.numDays);

Nx=(endLon-startLon)/dLon;
Ny=(endLat-startLat)/dLat;

Nt=bin.numDaysTot*NtPerDay;
bin.dt=bin.numDaysTot*24*60*60/Nt;

% Load MJO indices and read in ascat start date for lookup.
MJOdata=csvread('MJOdata.txt');
ascatStart=[1990 01 01 0 0 0];
   
% Adjust so average of bin is center of bin.
bin.x=startLon+dLon/2:dLon:endLon-dLon/2;
bin.y=startLat+dLat/2:dLat:endLat-dLat/2;

if(~(Nx == floor(Nx)) || ~(Ny == floor(Ny)) || ~(bin.numDaysTot == floor(bin.numDaysTot)))
    error('Invalid ranges or increments');
end

% u, v and time will contain bin averages
bin.p=zeros(length(bin.x),length(bin.y),Nt);
bin.time=zeros(length(bin.x),length(bin.y),Nt);

bin.count=zeros(length(bin.x),length(bin.y),Nt);

% MJOstring=num2str(binArgs.mjoPhases(1));
% for i=2:length(binArgs.mjoPhases)
%     MJOstring=[MJOstring ', ' num2str(binArgs.mjoPhases(i))];
% end

%--------------------------------------------------------------------------
% Bin the data. 
%--------------------------------------------------------------------------    

fprintf('Binning.\n');

% Loop through time intervals
for i=1:numIntervals
    
    fprintf('Interval %i of %i. \n',i,numIntervals);

    netCDF=struct;
    startTime=dateCell{2*i-1};
    
    staticFileName=findFiles('static.nc');
    
    % u, v and time will contain bin averages, I indicates interval i.
    pI=zeros(length(bin.x),length(bin.y),bin.numDays(i)*NtPerDay);
    countI=zeros(length(bin.x),length(bin.y),bin.numDays(i)*NtPerDay);
    timeI=zeros(length(bin.x),length(bin.y),bin.numDays(i)*NtPerDay);

    for k=1:length(fileNames{i})

        netCDF.filename=fileNames{i}{k};
        
        fprintf('File %i of %i. \n',k,length(fileNames{i}));
        
        % Obtain date from filename
        year=str2double(fileNames{i}{k}(length(fileNames{i}{k})-10:length(fileNames{i}{k})-7));
        month=str2double(fileNames{i}{k}(length(fileNames{i}{k})-6:length(fileNames{i}{k})-5));
        day=str2double(fileNames{i}{k}(length(fileNames{i}{k})-4:length(fileNames{i}{k})-3));
        netCDF.date=[year, month, day, 0 0 0];
        
        netCDF.lat=repmat(ncread(staticFileName{1},'XLAT',[1 1 1], [Inf Inf 1], [1 1 1]),[1,1,24]);
        netCDF.lon=repmat(ncread(staticFileName{1},'XLONG',[1 1 1], [Inf Inf 1], [1 1 1]),[1,1,24]);
        netCDF.p=squeeze(ncread(fileNames{i}{k},'RAINNC',[1 1 1], [Inf Inf Inf], [1 1 1]))*100*10;

        % These temp arrays will store the averaged lat and lon data.
        pSumLat=zeros(size(netCDF.lat,1),length(bin.y),24);
        tSumLat=zeros(size(netCDF.lat,1),length(bin.y),24);
        countLat=zeros(size(netCDF.lat,1),length(bin.y),24);

        % Form time array
        netCDF.time=ones(size(netCDF.p,1),size(netCDF.p,2),24);
        for n=1:24
            netCDF.time(:,:,n)=(n-1)*ones(size(netCDF.p,1),size(netCDF.p,2))*60*60;
        end

        % Transform to LST if required. 
        if binArgs.LST==1
            netCDF.time=netCDF.time+(netCDF.lon./360)*24*60*60;
        end

        % Carry out binning and averaging.
        kLat=floor((netCDF.lat(1,:,1)-startLat)/dLat)+1;
        kLon=floor((netCDF.lon(:,1,1)-startLon)/dLon)+1;
        kLon=repmat(kLon,[1,24]);
        kT=floor((etime(netCDF.date,startTime)+squeeze(netCDF.time(:,1,:)))/bin.dt)+1;

        % Time since startTime
        time=mod((netCDF.time(:,:,:))/3600-(startTime(4)+startTime(5)/60+startTime(6)/3600),24);
        
        % Sum into latitude bins first. 
        for m=1:length(bin.y)
            pSumLat(:,m,:)=sum(netCDF.p(:,kLat==m,:,:),2);
            tSumLat(:,m,:)=sum(time(:,kLat==m,:,:),2);
            countLat(:,m,:)=sum(kLat==m);
        end

        % Iterate through longitude and time bins.
            
        % Just work out the indices for the data we want to bin.
        [iInd, kInd] = find(kLon>=1 & kLon <= length(bin.x) & ...
            kT>=1 & kT<=bin.numDays(i)*NtPerDay);

        for m=1:length(iInd) % Iterate over swath and time.

            pI(kLon(iInd(m),kInd(m)),:,kT(iInd(m),kInd(m)))=...
                pI(kLon(iInd(m),kInd(m)),:,kT(iInd(m),kInd(m)))+pSumLat(iInd(m),:,kInd(m));
            timeI(kLon(iInd(m),kInd(m)),:,kT(iInd(m),kInd(m)))=...
                timeI(kLon(iInd(m),kInd(m)),:,kT(iInd(m),kInd(m)))+tSumLat(iInd(m),:,kInd(m));
            countI(kLon(iInd(m),kInd(m)),:,kT(iInd(m),kInd(m)))=...
                countI(kLon(iInd(m),kInd(m)),:,kT(iInd(m),kInd(m)))+countLat(iInd(m),:,kInd(m));
            
        end
    end
    
    % Add data from interval I to total u,v,count arrays.
    if i==1
        bin.p(:,:,1:bin.numDays(i)*NtPerDay)=pI;
        bin.count(:,:,1:bin.numDays(i)*NtPerDay)=countI;
        bin.time(:,:,1:bin.numDays(i)*NtPerDay)=timeI;

    else
        bin.p(:,:,sum(bin.numDays(1:i-1))*NtPerDay+1:sum(bin.numDays(1:i))*NtPerDay)=pI;
        bin.count(:,:,sum(bin.numDays(1:i-1))*NtPerDay+1:sum(bin.numDays(1:i))*NtPerDay)=countI;
        bin.time(:,:,sum(bin.numDays(1:i-1))*NtPerDay+1:sum(bin.numDays(1:i))*NtPerDay)=timeI;
    end
    
end

% If using the ASCAT method, only take values corresponding to the
% satellite track.
if ~isnan(ASCATcount)
    bin.p(ASCATcount==1)=0;
    bin.count(ASCATcount==1)=0;
    bin.time(ASCATcount==1)=0;
end

bin.p=bin.p./bin.count;
bin.time=bin.time./bin.count;

bin.p(isnan(bin.p) | isinf(bin.p))=0;
bin.time(isnan(bin.time) | isinf(bin.time))=0;

clear uBinI vBinI countI timeBinI mjoTest netCDF uSumLon uSumLonLat ...
    vSumLon vSumLonLat countLon countLonLat;

% Create coherent data structure for output
bin.dateCell=dateCell;

bin.NtPerDay=NtPerDay;
bin.binArgs=binArgs;

end

