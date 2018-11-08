function pValue = rankSumRSCAT(bin)

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
    
    data1=squeeze(bin.p(iInd(i),jInd(i),1:6:end));
    countTest=squeeze(bin.count(iInd(i),jInd(i),1:6:end));
    data1(countTest==0)=[];
    
    data2=squeeze(bin.p(iInd(i),jInd(i),2:6:end));
    countTest=squeeze(bin.count(iInd(i),jInd(i),2:6:end));
    data2(countTest==0)=[];
    
    data3=squeeze(bin.p(iInd(i),jInd(i),3:6:end));
    countTest=squeeze(bin.count(iInd(i),jInd(i),3:6:end));
    data3(countTest==0)=[];
    
    data4=squeeze(bin.p(iInd(i),jInd(i),4:6:end));
    countTest=squeeze(bin.count(iInd(i),jInd(i),4:6:end));
    data4(countTest==0)=[];
    
    data5=squeeze(bin.p(iInd(i),jInd(i),5:6:end));
    countTest=squeeze(bin.count(iInd(i),jInd(i),5:6:end));
    data5(countTest==0)=[];
    
    data6=squeeze(bin.p(iInd(i),jInd(i),6:6:end));
    countTest=squeeze(bin.count(iInd(i),jInd(i),6:6:end));
    data6(countTest==0)=[];
    
    try
        p1=ranksum(data1,data4);
        p2=ranksum(data2,data5);
        p3=ranksum(data3,data6);
        
        p4=ranksum(data1,data2);
        p5=ranksum(data1,data3);
        p6=ranksum(data1,data5);
        p7=ranksum(data1,data6);
        
        p8=ranksum(data2,data3);
        p9=ranksum(data2,data4);
        p10=ranksum(data2,data6);
        
        p11=ranksum(data3,data4);
        p12=ranksum(data3,data5);
        
        p13=ranksum(data4,data5);
        p14=ranksum(data4,data6);
        
        p15=ranksum(data5,data6);
        
        p=min([p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15]);
        pValue(iInd(i),jInd(i))=p;
    catch
        pValue(iInd(i),jInd(i))=nan;
        failureCount=failureCount+1;
    end
    
end

end

