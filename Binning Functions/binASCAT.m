function bin=binASCAT(fileNames,dLat,startLat,endLat,...
    dLon,startLon,endLon,NtPerDay,dateCell,binArgs)
% binASCAT.m ----------------------------------------------------------
% 
% Copyright Ewan Short 23/2/2016
%
% DESC: This function takes wind direction from netCDF format at time point in 
% column t and plots.
%
% INPUTS: fileNames is an array containing the fileNames of the netCDF files, 
% seperated into seperate intervals by row.
% dLat and dLon are the
% size of lat and lon increments respectively, dt time increment in 
% seconds. Note the
% time variable in the netCDF data is measured in seconds since 1/1/1990
% midnight. tRange is a vector [startTime endTime] with components in
% same units. If averageTimes is 1, average data for given time bins across
% all days in data range. 
%
% OUTPUT: Output the binned average u and v component winds and a count of 
% how many measurements were in each bin for averaging purposes. 
%
%--------------------------------------------------------------------------

if nargin<10
    binArgs=struct;
    % Use LST (1) or UTC (0).
    binArgs.LST=1; 
    % Only bin certain MJO phases (1)
    binArgs.mjo=0;
    % Which phases?
    binArgs.mjoPhases=[1 2 3 4 5 6 7 8];
    % Set number of WRF data points to skip.
    binArgs.dataSource='ASCAT 12.5km Coast Opt.';
end

close all;

% Determine number of time intervals. Check have been inputted correctly.
numTimes=length(dateCell);
numIntervals=numTimes/2;

if floor(numIntervals)~=numIntervals
    error('Odd number of date vectors. \n');
end

if length(binArgs.mjoPhases)<8 && binArgs.mjo==0
    warning('binArgs.mjo==0');
end
    
% Calculate total number of days and confirm an integer. 
bin.numDays=zeros(1,numIntervals);

for i=1:(numIntervals)
    bin.numDays(i)=etime(dateCell{2*i},dateCell{2*i-1})/(24*60*60);
    if floor(bin.numDays(i))~=bin.numDays
        error('Time intervals must be in whole days. \n');
    end
end

% Confirm number of intervals same as number of rows of fileNames
if (size(fileNames,1)~=numIntervals)
    error('Size of fileNames cell does not match number of intervals');
end

bin.numDaysTot=sum(bin.numDays);

Nx=(endLon-startLon)/dLon;
Ny=(endLat-startLat)/dLat;

Nt=bin.numDaysTot*NtPerDay;
bin.dt=bin.numDaysTot*24*60*60/Nt;
   
% Adjust so average of bin is center of bin.
bin.x=startLon+dLon/2:dLon:endLon-dLon/2;
bin.y=startLat+dLat/2:dLat:endLat-dLat/2;
% Recall ascat time values are measured in seconds since 1/1/1990 midnight UTC.
% Include option to add leap seconds manually. 
ascatStart=[1990 01 01 0 0 0];
% Read in MJO indices
MJOdata=csvread('MJOdata.txt');

if(~(Nx == floor(Nx)) || ~(Ny == floor(Ny)))
    error('Invalid ranges or increments');
end

bin.u=zeros(length(bin.x),length(bin.y),Nt);
bin.v=zeros(length(bin.x),length(bin.y),Nt);
bin.time=zeros(length(bin.x),length(bin.y),Nt);
bin.count=zeros(length(bin.x),length(bin.y),Nt);

%--------------------------------------------------------------------------
% Bin the data. 
%--------------------------------------------------------------------------        
    
fprintf('Binning.\n');

% Loop through time intervals
for i=1:numIntervals
    
    fprintf('Interval %i of %i. \n',i,numIntervals);

    netCDF=struct;
    startTime=dateCell{2*i-1};
    endTime=dateCell{2*i};
    
    uI=zeros(length(bin.x),length(bin.y),bin.numDays(i)*NtPerDay);
    vI=zeros(length(bin.x),length(bin.y),bin.numDays(i)*NtPerDay);
    timeI=zeros(length(bin.x),length(bin.y),bin.numDays(i)*NtPerDay);
    countI=zeros(length(bin.x),length(bin.y),bin.numDays(i)*NtPerDay,'uint8');

    for k=1:length(fileNames{i})

        netCDF.filename=fileNames{i}{k};
        netCDF.lat=ncread(fileNames{i}{k},'lat');
        netCDF.lon=ncread(fileNames{i}{k},'lon');
        netCDF.time=ncread(fileNames{i}{k},'time');
        netCDF.wind_speed=ncread(fileNames{i}{k},'wind_speed');
        netCDF.wind_dir=ncread(fileNames{i}{k},'wind_dir');

        year=str2double(fileNames{i}{k}(length(fileNames{i}{k})-52:length(fileNames{i}{k})-49));
        month=str2double(fileNames{i}{k}(length(fileNames{i}{k})-48:length(fileNames{i}{k})-47));
        day=str2double(fileNames{i}{k}(length(fileNames{i}{k})-46:length(fileNames{i}{k})-45));
        netCDF.date=[year, month, day, 0 0 0];

        if binArgs.mjo
            % Specify MJO index from BOM data
            netCDF.mjoPhase=MJOdata(etime(netCDF.date,ascatStart)./(24*60*60)+1,4);
            netCDF.mjoAmp=MJOdata(etime(netCDF.date,ascatStart)./(24*60*60)+1,5);

            % Check if correct MJO phase and amplitude
            mjoTest=(any(netCDF.mjoPhase==binArgs.mjoPhases) && netCDF.mjoAmp>=1)*...
                ones(size(netCDF.lat,1),size(netCDF.lat,2));
        else
            mjoTest=ones(size(netCDF.lat,1),size(netCDF.lat,2));
        end

        % If using local solar time, adjust seconds since midnight 1/1/1990 UTC
        % accordingly. Does it make sense to add this number? When should be
        % subtract it? Below formula makes sense when lon expressed in +/-
        % 180 convection. Need to change longitudes >180... i.e subtract
        % 360 from them...?
        if binArgs.LST==1
            netCDF.time=netCDF.time+(netCDF.lon/360)*24*60*60;
        end

        % Calculate time constraint.
        startTimeSec=etime(startTime,ascatStart);
        endTimeSec=etime(endTime,ascatStart);
        startTimeTest=netCDF.time-startTimeSec;
        endTimeTest=endTimeSec-netCDF.time;

        kLat=floor((netCDF.lat-startLat)/dLat)+1;
        kLon=floor((netCDF.lon-startLon)/dLon)+1;
        kT=floor((etime(ascatStart,startTime)+netCDF.time)/bin.dt)+1;
        
        % Just work out the indices for the data we want to bin.
        [iInd, jInd] = find(kLat >= 1 & kLat <= length(bin.y) ...
            & kLon>=1 & kLon <= length(bin.x) & ...
            ~isnan(netCDF.wind_speed) & ~isnan(netCDF.wind_dir) & ...
            kT>=1 & kT<=bin.numDays(i)*NtPerDay & mjoTest>0);
        
        for m=1:length(iInd) % Iterate over swath and time.

            uI(kLon(iInd(m),jInd(m)),kLat(iInd(m),jInd(m)),kT(iInd(m),jInd(m)))=...
                uI(kLon(iInd(m),jInd(m)),kLat(iInd(m),jInd(m)),kT(iInd(m),jInd(m)))+netCDF.wind_speed(iInd(m),jInd(m))*...
                sind(netCDF.wind_dir(iInd(m),jInd(m)));
            vI(kLon(iInd(m),jInd(m)),kLat(iInd(m),jInd(m)),kT(iInd(m),jInd(m)))=...
                vI(kLon(iInd(m),jInd(m)),kLat(iInd(m),jInd(m)),kT(iInd(m),jInd(m)))+netCDF.wind_speed(iInd(m),jInd(m))*...
                cosd(netCDF.wind_dir(iInd(m),jInd(m)));
            timeI(kLon(iInd(m),jInd(m)),kLat(iInd(m),jInd(m)),kT(iInd(m),jInd(m)))=...
                timeI(kLon(iInd(m),jInd(m)),kLat(iInd(m),jInd(m)),kT(iInd(m),jInd(m)))...
                +rem(startTimeTest(iInd(m),jInd(m)),24*60*60)/(60*60);

            % Add one to number of measurements in this bin. 
            countI(kLon(iInd(m),jInd(m)),kLat(iInd(m),jInd(m)),kT(iInd(m),jInd(m)))=...
                countI(kLon(iInd(m),jInd(m)),kLat(iInd(m),jInd(m)),kT(iInd(m),jInd(m)))+1;        
        end

    end
    
    if i==1
        bin.u(:,:,1:bin.numDays(i)*NtPerDay)=uI;
        bin.v(:,:,1:bin.numDays(i)*NtPerDay)=vI;
        bin.count(:,:,1:bin.numDays(i)*NtPerDay)=countI;
        bin.time(:,:,1:bin.numDays(i)*NtPerDay)=timeI;
    else
        bin.u(:,:,sum(bin.numDays(1:i-1))*NtPerDay+1:sum(bin.numDays(1:i))*NtPerDay)=uI;
        bin.v(:,:,sum(bin.numDays(1:i-1))*NtPerDay+1:sum(bin.numDays(1:i))*NtPerDay)=vI;
        bin.count(:,:,sum(bin.numDays(1:i-1))*NtPerDay+1:sum(bin.numDays(1:i))*NtPerDay)=countI;
        bin.time(:,:,sum(bin.numDays(1:i-1))*NtPerDay+1:sum(bin.numDays(1:i))*NtPerDay)=timeI;
    end
    
    clear uBinI vBinI countI timeBinI netCDF startTimeSec endTimeSec startTimeTest...
        endTimeTest mjoTest year month day mjoData
    
end

%--------------------------------------------------------------------------
% Average observations in each bin
%--------------------------------------------------------------------------

bin.u(bin.count>0)=bin.u(bin.count>0)./bin.count(bin.count>0);
bin.v(bin.count>0)=bin.v(bin.count>0)./bin.count(bin.count>0);
bin.time(bin.count>0)=bin.time(bin.count>0)./bin.count(bin.count>0);

% Create coherent data structure for output
bin.dateCell=dateCell;
bin.NtPerDay=NtPerDay;
bin.binArgs=binArgs;

end

