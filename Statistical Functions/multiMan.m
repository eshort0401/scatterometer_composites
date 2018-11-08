function pValue = multiMan(bin)

failureCount=0;

pValue=nan(length(bin.x),length(bin.y));

% Check if countComp exists and if not create it
if ~isfield(bin,'countComp')
    % Create array to store composite count.    
    bin.countComp=zeros(length(bin.x),length(bin.y),bin.NtPerDay);

    for i=1:bin.NtPerDay
        bin.countComp(:,:,i)=sum(bin.count(:,:,i+(0:bin.numDaysTot-1)*bin.NtPerDay)>0,3);
    end 
end

% Find spatial points with data
[iInd, jInd] = find(all(bin.countComp,3)>0 & all(bin.countComp,3)>0);

for i=1:length(iInd)
    
    % Now were only working with the spatial point (iInd(i),jInd(i)).
    % First create array of the form (u,v,t) for this point.
    
    kInd=find(bin.uPert(iInd(i),jInd(i),:)~=0);
    kIndMod=mod(kInd-1,12)+1;
    
    X=[squeeze(bin.uPert(iInd(i),jInd(i),kInd)) squeeze(bin.vPert(iInd(i),jInd(i),kInd)) kIndMod];
    [~, I]=sort(X(:,3));
    X=X(I,:);
    
    try
        [~, p]=manova1(X(:,1:2),X(:,3));
        pValue(iInd(i),jInd(i))=p(1);
    catch
        pValue(iInd(i),jInd(i))=nan;
        failureCount=failureCount+1;
    end
    
end

end

