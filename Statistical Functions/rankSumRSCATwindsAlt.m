function pValue = rankSumRSCATwindsAlt(bin)

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
      
    data2u=squeeze(bin.u(iInd(i),jInd(i),2:6:end));
    countTest=squeeze(bin.count(iInd(i),jInd(i),2:6:end));
    data2u(countTest==0)=[];
    
    data2v=squeeze(bin.v(iInd(i),jInd(i),2:6:end));
    data2v(countTest==0)=[];
    
    data5u=squeeze(bin.u(iInd(i),jInd(i),5:6:end));
    countTest=squeeze(bin.count(iInd(i),jInd(i),5:6:end));
    data5u(countTest==0)=[];
    
    data5v=squeeze(bin.v(iInd(i),jInd(i),5:6:end));
    data5v(countTest==0)=[];
    
    try
        p2u=ranksum(data2u,data5u);
        p2v=ranksum(data2v,data5v);
                
        p=min([p2u,p2v]);
        pValue(iInd(i),jInd(i))=p;
    catch
        pValue(iInd(i),jInd(i))=nan;
        failureCount=failureCount+1;
    end
    
end

end

