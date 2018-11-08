function [pValue] = rankSumASCAT(bin1,bin2,T1,T2)

failureCount=0;

pValue=nan(length(bin1.x),length(bin1.y));
pa1i1=nan(length(bin1.x),length(bin1.y));
pa1a2=nan(length(bin1.x),length(bin1.y));
pa2i2=nan(length(bin1.x),length(bin1.y));
pi1i2=nan(length(bin1.x),length(bin1.y));

% Check if countComp exists and if not create it
if ~isfield(bin1,'countComp')
    % Create array to store composite count.    
    bin1.countComp=zeros(length(bin1.x),length(bin1.y),bin1.NtPerDay);

    for i=1:bin1.NtPerDay
        bin1.countComp(:,:,i)=sum(bin1.count(:,:,i+(0:bin1.numDaysTot-1)*bin1.NtPerDay)>0,3);
        bin2.countComp(:,:,i)=sum(bin2.count(:,:,i+(0:bin1.numDaysTot-1)*bin2.NtPerDay)>0,3);
    end 
end

% Find spatial points with data
[iInd, jInd] = find(bin1.countComp(:,:,T1)>0 & bin1.countComp(:,:,T2)>0 & ...
    bin2.countComp(:,:,T1)>0 & bin2.countComp(:,:,T2)>0);

for i=1:length(iInd)
    
    % Now were only working with the spatial point (iInd(i),jInd(i)).
    % First create array of the form (u,v,t) for this point.
    
    data1=squeeze(bin1.p(iInd(i),jInd(i),T1:16:end));
    countTest=squeeze(bin1.count(iInd(i),jInd(i),T1:16:end));
    data1(countTest==0)=[];
    
    data2=squeeze(bin2.p(iInd(i),jInd(i),T1:16:end));
    countTest=squeeze(bin2.count(iInd(i),jInd(i),T1:16:end));
    data2(countTest==0)=[];
    
    data3=squeeze(bin1.p(iInd(i),jInd(i),T2:16:end));
    countTest=squeeze(bin1.count(iInd(i),jInd(i),T2:16:end));
    data3(countTest==0)=[];
    
    data4=squeeze(bin2.p(iInd(i),jInd(i),T2:16:end));
    countTest=squeeze(bin2.count(iInd(i),jInd(i),T2:16:end));
    data4(countTest==0)=[];

    try
        p1=ranksum(data1,data2); % Active T1 with inactive T1
        p2=ranksum(data1,data3); % Active T1 with Active T2
        p3=ranksum(data2,data4); % Inactive T1 with Inactive T2
        p4=ranksum(data3,data4); % Active T2 with Inactive T2
        
        pa1i1(iInd(i),jInd(i))=p1; % Active T1 with inactive T1
        pa1a2(iInd(i),jInd(i))=p2; % Active T1 with Active T2
        pi1i2(iInd(i),jInd(i))=p3; % Inactive T1 with Inactive T2
        pa2i2(iInd(i),jInd(i))=p4; % Active T2 with Inactive T2
        
        p=min([p1,p4]);
        pValue(iInd(i),jInd(i))=p;
    catch
        pValue(iInd(i),jInd(i))=nan;
        failureCount=failureCount+1;
    end
    
end

end

