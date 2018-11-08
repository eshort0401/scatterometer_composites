function bin=binHY2SCAT(fileNames,dLat,startLat,endLat,...
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
    binArgs.dataSource='HY2SCAT 25km.';
    binArgs.corruptedFiles={};
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

    HDF5=struct;
    startTime=dateCell{2*i-1};
    endTime=dateCell{2*i};
    
    uI=zeros(length(bin.x),length(bin.y),bin.numDays(i)*NtPerDay);
    vI=zeros(length(bin.x),length(bin.y),bin.numDays(i)*NtPerDay);
    timeI=zeros(length(bin.x),length(bin.y),bin.numDays(i)*NtPerDay);
    countI=zeros(length(bin.x),length(bin.y),bin.numDays(i)*NtPerDay,'uint8');
    
    corruptCount=1;

    for k=1:length(fileNames{i})

        try
            HDF5.filename=fileNames{i}{k};
            % Note the HY2SCAT arrays only provide lat,lon values where there
            % is data.
            HDF5.lat=h5read(fileNames{i}{k},'/wvc_lat');
            HDF5.lat(HDF5.lat==0)=NaN;
            HDF5.lon=h5read(fileNames{i}{k},'/wvc_lon');
            HDF5.lon(HDF5.lon==0)=NaN;
            HDF5.time=h5read(fileNames{i}{k},'/row_time');
        catch
            warning('File Corrupted.');
            binArgs.corruptedFiles{corruptCount}=fileNames{i}{k};
            corruptCount=corruptCount+1;
            continue;
        end
            
        % Look for nonzero values
        
        
        for j=1:1624
            if HDF5.time{j}(9)=='T'
                HDF5.time{j}=[str2double(HDF5.time{j}(1:4)) ...
                    str2double(HDF5.time{j}(5:6)) str2double(HDF5.time{j}(7:8)) ...
                    str2double(HDF5.time{j}(10:11)) str2double(HDF5.time{j}(13:14)) ...
                    str2double(HDF5.time{j}(16:17))];
                    HDF5.time{j}=etime(HDF5.time{j},refTime);
            else
                HDF5.time{j}=NaN;
            end
        end
        HDF5.time=cell2mat(HDF5.time);
        HDF5.time=repmat(HDF5.time',76,1);
        
        HDF5.wind_speed=double(h5read(fileNames{i}{k},'/wind_speed_selection'))*0.01;
        HDF5.wind_dir=double(h5read(fileNames{i}{k},'/wind_dir_selection'))*0.1;

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
        kLon=floor((HDF5.lon-startLon)/dLon)+1;
        kT=floor((etime(refTime,startTime)+HDF5.time)/bin.dt)+1;
        
        % Just work out the indices for the data we want to bin.
        [iInd, jInd] = find(kLat >= 1 & kLat <= length(bin.y) ...
            & kLon>=1 & kLon <= length(bin.x) & ...
            ~isnan(HDF5.wind_speed) & ~isnan(HDF5.wind_dir) & ...
            kT>=1 & kT<=bin.numDays(i)*NtPerDay);
        
        for m=1:length(iInd) % Iterate over swath and time.

            uI(kLon(iInd(m),jInd(m)),kLat(iInd(m),jInd(m)),kT(iInd(m),jInd(m)))=...
                uI(kLon(iInd(m),jInd(m)),kLat(iInd(m),jInd(m)),kT(iInd(m),jInd(m)))+HDF5.wind_speed(iInd(m),jInd(m))*...
                sind(HDF5.wind_dir(iInd(m),jInd(m)));
            vI(kLon(iInd(m),jInd(m)),kLat(iInd(m),jInd(m)),kT(iInd(m),jInd(m)))=...
                vI(kLon(iInd(m),jInd(m)),kLat(iInd(m),jInd(m)),kT(iInd(m),jInd(m)))+HDF5.wind_speed(iInd(m),jInd(m))*...
                cosd(HDF5.wind_dir(iInd(m),jInd(m)));
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
    
    clear uBinI vBinI countI timeBinI HDF5 startTimeSec endTimeSec startTimeTest...
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

