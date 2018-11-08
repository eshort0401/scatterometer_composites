function bin=binGPM(fileNames,dLat,startLat,endLat,...
    dLon,startLon,endLon,NtPerDay,dateCell,binArgs)
% plot_wind_dir.m ----------------------------------------------------------
% 
% Copyright Ewan Short 23/2/2016
%
% DESC: This function takes wind direction from HDF5 format at time point in 
% column t and plots.
%
% INPUTS: fileNames is an array containing the fileNames of the HDF5 files, 
% seperated into seperate intervals by row.
% dLat and dLon are the
% size of lat and lon increments respectively, dt time increment in 
% seconds. Note the
% time variable in the HDF5 data is measured in seconds since 1/1/1990
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
    % Set number of WRF data points to skip.
    binArgs.dataSource='GPM 25km.';
    binArgs.corruptedFiles={};
    binArgs.mjo=0;
    binArgs.mjoPhases=[];
    binArgs.mjoAmp=0;
end

close all;

% dLat=0.25;
% dLon=0.25;

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
% Use this as a reference time
refTime=[1990 01 01 0 0 0];

if(~(Nx == floor(Nx)) || ~(Ny == floor(Ny)))
    error('Invalid ranges or increments');
end

bin.p=zeros(length(bin.x),length(bin.y),Nt);
bin.time=zeros(length(bin.x),length(bin.y),Nt);
bin.count=zeros(length(bin.x),length(bin.y),Nt);

MJOdata=csvread('MJOdata.txt');

%--------------------------------------------------------------------------
% Bin the data. 
%--------------------------------------------------------------------------        
    
fprintf('Binning.\n');

% Loop through time intervals
for i=1:numIntervals
    
    fprintf('Interval %i of %i. \n',i,numIntervals);

    HDF5=struct;
    startTime=dateCell{2*i-1};
    endTime=dateCell{2*i};
    
    pI=zeros(length(bin.x),length(bin.y),bin.numDays(i)*NtPerDay);
    timeI=zeros(length(bin.x),length(bin.y),bin.numDays(i)*NtPerDay);
    countI=zeros(length(bin.x),length(bin.y),bin.numDays(i)*NtPerDay,'uint8');
    
    corruptCount=1;

    for k=1:length(fileNames{i})

        try
            HDF5.filename=fileNames{i}{k};
            % Note the HY2SCAT arrays only provide lat,lon values where there
            % is data.
            HDF5.lat=h5read(fileNames{i}{k},'/NS/Latitude');
            HDF5.lat(HDF5.lat==-9999.9)=NaN;
            HDF5.lon=h5read(fileNames{i}{k},'/NS/Longitude');
            HDF5.lon(HDF5.lon==-9999.9)=NaN;
            HDF5.year=double(h5read(fileNames{i}{k},'/NS/ScanTime/Year'));
            HDF5.month=double(h5read(fileNames{i}{k},'/NS/ScanTime/Month'));
            HDF5.day=double(h5read(fileNames{i}{k},'/NS/ScanTime/DayOfMonth'));
            HDF5.secondOfDay=double(h5read(fileNames{i}{k},'/NS/ScanTime/SecondOfDay'));
            timeN=length(HDF5.secondOfDay);
            HDF5.time=cell(1,timeN);
            HDF5.precip=double(h5read(fileNames{i}{k},'/NS/SLV/precipRateESurface'));
            % Set Latitude values to NaN wherever precip is corrupted for
            % easy rejection.
            HDF5.lat(HDF5.precip<0)=NaN;
            HDF5.mjoPhase=-99*ones(size(HDF5.time(1,:)));
            HDF5.mjoAmp=-99*ones(size(HDF5.time(1,:)));
            if timeN==0
                error('File Empty!');
            end

        catch
            warning('File Corrupted.');
            binArgs.corruptedFiles{corruptCount}=fileNames{i}{k};
            corruptCount=corruptCount+1;
            continue;
        end
                    
        for j=1:timeN
            if HDF5.secondOfDay(j)>=0
                HDF5.time{j}=[HDF5.year(j) HDF5.month(j) HDF5.day(j) 0 0 0];
                    HDF5.time{j}=etime(HDF5.time{j},refTime)+HDF5.secondOfDay(j);
            else
                HDF5.time{j}=NaN;
            end
        end
        HDF5.time=cell2mat(HDF5.time);
        HDF5.time=repmat(HDF5.time,49,1);
        
        if binArgs.mjo

            notNanIndex=find(~isnan(HDF5.time(1,:)));

            % Specify MJO index from BOM data
            HDF5.mjoPhase(notNanIndex)=MJOdata(floor(HDF5.time(1,notNanIndex)/(24*60*60))+1,4);
            HDF5.mjoAmp(notNanIndex)=MJOdata(floor(HDF5.time(1,notNanIndex)/(24*60*60))+1,5);

            % Check if correct MJO phase and amplitude
            mjoTest=(any(HDF5.mjoPhase'==binArgs.mjoPhases,2) & (HDF5.mjoAmp>=1)');

        else
            mjoTest=ones(timeN,1);
        end
        
        mjoTest=repmat(mjoTest,1,49)';
        
        % If using local solar time, adjust seconds since midnight 1/1/1990 UTC
        % accordingly. Does it make sense to add this number? When should be
        % subtract it?
        if binArgs.LST==1
            HDF5.time=HDF5.time+(HDF5.lon/360)*24*60*60;
        end

        % Calculate time constraint.
        startTimeSec=etime(startTime,refTime);
        endTimeSec=etime(endTime,refTime);
        startTimeTest=HDF5.time-startTimeSec;
        endTimeTest=endTimeSec-HDF5.time;

        kLat=floor((HDF5.lat-startLat)/dLat)+1;
        kLon=floor(mod(HDF5.lon-startLon,360)/dLon)+1;
        kT=floor((etime(refTime,startTime)+HDF5.time)/bin.dt)+1;
        
        % Just work out the indices for the data we want to bin.
        [iInd, jInd] = find(kLat >= 1 & kLat <= length(bin.y) ...
            & kLon>=1 & kLon <= length(bin.x) & ...
            ~isnan(HDF5.precip) & ...
            kT>=1 & kT<=bin.numDays(i)*NtPerDay & mjoTest>0);
        
        for m=1:length(iInd) % Iterate over swath and time.

            pI(kLon(iInd(m),jInd(m)),kLat(iInd(m),jInd(m)),kT(iInd(m),jInd(m)))=...
                pI(kLon(iInd(m),jInd(m)),kLat(iInd(m),jInd(m)),kT(iInd(m),jInd(m)))+HDF5.precip(iInd(m),jInd(m));
            % The time variable will store time in seconds since startTime
            % mod 24*60*60. This simplifies averaging later.
            timeI(kLon(iInd(m),jInd(m)),kLat(iInd(m),jInd(m)),kT(iInd(m),jInd(m)))=...
                timeI(kLon(iInd(m),jInd(m)),kLat(iInd(m),jInd(m)),kT(iInd(m),jInd(m)))...
                +rem(startTimeTest(iInd(m),jInd(m)),24*60*60)/(60*60);

            % Add one to number of measurements in this bin. 
            countI(kLon(iInd(m),jInd(m)),kLat(iInd(m),jInd(m)),kT(iInd(m),jInd(m)))=...
                countI(kLon(iInd(m),jInd(m)),kLat(iInd(m),jInd(m)),kT(iInd(m),jInd(m)))+1;        
        end

    end
    
    if i==1
        bin.p(:,:,1:bin.numDays(i)*NtPerDay)=pI;
        bin.count(:,:,1:bin.numDays(i)*NtPerDay)=countI;
        bin.time(:,:,1:bin.numDays(i)*NtPerDay)=timeI;
    else
        bin.p(:,:,sum(bin.numDays(1:i-1))*NtPerDay+1:sum(bin.numDays(1:i))*NtPerDay)=pI;
        bin.count(:,:,sum(bin.numDays(1:i-1))*NtPerDay+1:sum(bin.numDays(1:i))*NtPerDay)=countI;
        bin.time(:,:,sum(bin.numDays(1:i-1))*NtPerDay+1:sum(bin.numDays(1:i))*NtPerDay)=timeI;
    end
    
    clear pI countI timeBinI HDF5 startTimeSec endTimeSec startTimeTest...
        endTimeTest mjoTest year month day mjoData
    
end

%--------------------------------------------------------------------------
% Average observations in each bin
%--------------------------------------------------------------------------

bin.p(bin.count>0)=bin.p(bin.count>0)./bin.count(bin.count>0);
bin.time(bin.count>0)=bin.time(bin.count>0)./bin.count(bin.count>0);

% Create coherent data structure for output
bin.dateCell=dateCell;
bin.NtPerDay=NtPerDay;
bin.binArgs=binArgs;

end

