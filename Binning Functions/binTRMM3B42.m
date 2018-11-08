function bin = binTRMM3B42(fileNames,startLat,endLat,...
    startLon,endLon,NtPerDay,dateCell,binArgs)

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
    binArgs.dataSource='TRMM 3B42 25km';
    binArgs.corruptedFiles={};
end

close all;

dLat=0.25;
dLon=0.25;

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

%--------------------------------------------------------------------------
% Bin the data. 
%--------------------------------------------------------------------------        
    
fprintf('Binning.\n');

% Loop through time intervals
for i=1:numIntervals
    
    fprintf('Interval %i of %i. \n',i,numIntervals);

    HDF=struct;
    startTime=dateCell{2*i-1};
    endTime=dateCell{2*i};
    HDF.lon=(-180+.125):.25:(180-.125);
    HDF.lon=repmat(HDF.lon',[1 400]);
    HDF.lat=(-50+.125):.25:(50-.125);
    HDF.lat=repmat(HDF.lat,[1440 1]);
    
    pI=zeros(length(bin.x),length(bin.y),bin.numDays(i)*NtPerDay);
    timeI=zeros(length(bin.x),length(bin.y),bin.numDays(i)*NtPerDay);
    countI=zeros(length(bin.x),length(bin.y),bin.numDays(i)*NtPerDay,'uint8');
    
    corruptCount=1;

    for k=1:length(fileNames{i})

        try
            HDF.fileHeader=hdfread(fileNames{i}{k},'FileHeader');
            HDF.fileHeader=HDF.fileHeader{1}';
            % Note the HY2SCAT arrays only provide lat,lon values where there
            % is data.
            HDF.satObTime=double(hdfread(fileNames{i}{k},'satObservationTime'));
            % This is measured in mm/h.
            HDF.precip=hdfread(fileNames{i}{k},'precipitation');
            HDF.time=[str2double(HDF.fileHeader(60:63))...
                str2double(HDF.fileHeader(64:65))...
                str2double(HDF.fileHeader(66:67))...
                str2double(HDF.fileHeader(69:70)) 0 0];
            % TRMM data uses plus minus 180 convention
            
            % Use local solar time
            
            if binArgs.LST
                HDF.time=etime(HDF.time,refTime)+HDF.satObTime*60+...
                    (HDF.lon/360)*24*60*60;
            else
                HDF.time=etime(HDF.time,refTime)+HDF.satObTime*60;
            end
        catch
            warning('File Corrupted.');
            binArgs.corruptedFiles{corruptCount}=fileNames{i}{k};
            corruptCount=corruptCount+1;
            continue;
        end

        % Calculate time constraint.
        startTimeSec=etime(startTime,refTime);
        endTimeSec=etime(endTime,refTime);
        startTimeTest=HDF.time-startTimeSec;
        endTimeTest=endTimeSec-HDF.time;

        kLat=floor((HDF.lat-startLat)/dLat)+1;
        kLon=floor((HDF.lon-startLon)/dLon)+1;
        kT=floor((etime(refTime,startTime)+HDF.time)/bin.dt)+1;
        
        % Just work out the indices for the data we want to bin.
        [iInd, jInd] = find(kLat >= 1 & kLat <= length(bin.y) ...
            & kLon>=1 & kLon <= length(bin.x) & ...
            HDF.precip>0 & HDF.precip~=-9999.9 & ...
            kT>=1 & kT<=bin.numDays(i)*NtPerDay);
        
        for m=1:length(iInd) % Iterate over lat, lon

            pI(kLon(iInd(m),jInd(m)),kLat(iInd(m),jInd(m)),kT(iInd(m),jInd(m)))=...
                pI(kLon(iInd(m),jInd(m)),kLat(iInd(m),jInd(m)),kT(iInd(m),jInd(m)))...
                +HDF.precip(iInd(m),jInd(m));
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
        bin.p(:,:,sum(bin.numDays(1:i-1))*NtPerDay+1:sum(bin.numDays(1:i))*NtPerDay)=uI;
        bin.count(:,:,sum(bin.numDays(1:i-1))*NtPerDay+1:sum(bin.numDays(1:i))*NtPerDay)=countI;
        bin.time(:,:,sum(bin.numDays(1:i-1))*NtPerDay+1:sum(bin.numDays(1:i))*NtPerDay)=timeI;
    end
    
    clear pI countI timeI HDF startTimeSec endTimeSec startTimeTest...
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

