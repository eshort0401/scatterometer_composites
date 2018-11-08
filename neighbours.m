function [index, closest] = neighbours(count,NtPerDay,avRan,times,dateCell)

indexCount=1;

indexLength=sum(sum(sum(count>0)));
forward=zeros(indexLength,length(times),'uint16');
backward=zeros(indexLength,length(times),'uint16');
index=zeros(indexLength,3,'uint16');
kSubList=(1:size(count,3))';

kInd=find(any(any(count>0,1),2));

% Determine number of time intervals. Check have been inputted correctly.
numTimes=length(dateCell);
numIntervals=numTimes/2;
    
% Calculate total number of days and confirm an integer. 
numDays=zeros(1,numIntervals);
jump=zeros(1,numIntervals-1);

for i=1:numIntervals
    numDays(i)=etime(dateCell{2*i},dateCell{2*i-1})/(24*60*60);
    if floor(numDays(i))~=numDays
        error('Time intervals must be in whole days. \n');
    end
end

for i=1:(numIntervals-1)
    jump(i)=etime(dateCell{2*i+1},dateCell{2*i})/(24*60*60);
    if floor(numDays(i))~=numDays
        error('Time intervals must be in whole days. \n');
    end
end

for p=1:length(kInd)
    
    k=kInd(p);
    
    % Only look at nonzero bins
    [iInd, jInd]=find(count(:,:,k)>0);

    for m=1:length(iInd)
        
        i=iInd(m);
        j=jInd(m);
        
        index(indexCount,:)=[i j k];
        ijCount(:)=count(i,j,:);
                   
        % Start at current time, look for the next avRan "complete days"
        for n=1:length(times)
            
            l=times(n)-1;
            lCountForward=ijCount(k+l:NtPerDay:end);
            lCountBackward=ijCount(k-l:-NtPerDay:1);
            
            kSubListForward=kSubList(k+l:NtPerDay:end);
            kSubListBackward=kSubList(k-l:-NtPerDay:1);
            
            % Find the k index for the nearest avRan non-zero values
            sub=find(lCountForward>0,avRan,'first');
            kSubTemp=kSubListForward(sub);
            kSub=[kSubTemp zeros(1,avRan-length(kSubTemp))];
            forward(indexCount,n:length(times):length(times)*avRan)=kSub;
            
            sub=find(lCountBackward>0,avRan,'first');
            kSubTemp=kSubListBackward(sub);
            kSub=[kSubTemp zeros(1,avRan-length(kSubTemp))];
            backward(indexCount,n:length(times):length(times)*avRan)=kSub;
            
        end
      
        indexCount=indexCount+1;
        
    end
    
end

forward=single(forward);
backward=single(backward);
forward(forward==0)=nan;
backward(backward==0)=nan;
index=single(index);

% Add jump to forward, backward, index to account for disjoint intervals
% I.e the k values will actually reflect the distance in time between bins
forwardMod=zeros(indexLength,length(times));
forwardMod(forward<=numDays(1)*NtPerDay)=forward(forward<=numDays(1)*NtPerDay);
forwardMod(isnan(forward))=NaN;

backwardMod=zeros(indexLength,length(times));
backwardMod(backward<=numDays(1)*NtPerDay)=backward(backward<=numDays(1)*NtPerDay);
backwardMod(isnan(backward))=NaN;

indexMod=zeros(size(index));
indexMod(:,[1 2])=index(:,[1 2]);
indexMod(index(:,3)<=numDays(1)*NtPerDay)=index(index(:,3)<=numDays(1)*NtPerDay);
indexMod(isnan(index))=NaN;

for i=1:length(numDays)-1
    
    forwardIn=forward>sum(numDays(1:i))*NtPerDay & forward<=sum(numDays(1:i+1))*NtPerDay;
    backwardIn=backward>sum(numDays(1:i))*NtPerDay & backward<=sum(numDays(1:i+1))*NtPerDay;
    indexIn=index(:,3)>sum(numDays(1:i))*NtPerDay & index(:,3)<=sum(numDays(1:i+1))*NtPerDay;
    
    forwardMod(forwardIn)=forward(forwardIn)+sum(jump(1:i))*NtPerDay;
    backwardMod(backwardIn)=backward(backwardIn)+sum(jump(1:i))*NtPerDay;
    indexMod(indexIn,3)=index(indexIn,3)+sum(jump(1:i))*NtPerDay;

end

forDis=forwardMod-repmat(indexMod(:,3),1,length(times));
backDis=repmat(indexMod(:,3),1,length(times))-backwardMod;

closest=nan(indexLength,length(times),'single');
sansNanBoth=find(~isnan(forward) & ~isnan(backward));
sansNanFor=find(~isnan(forward) & isnan(backward));
sansNanBack=find(~isnan(backward) & isnan(forward));
closest(sansNanBoth)=forward(sansNanBoth)+...
    (backDis(sansNanBoth)<forDis(sansNanBoth)).*...
    (backward(sansNanBoth)-forward(sansNanBoth));
closest(sansNanFor)=forward(sansNanFor);
closest(sansNanBack)=backward(sansNanBack);

index(any(isnan(closest),2),:)=[];
closest(any(isnan(closest),2),:)=[];

end

