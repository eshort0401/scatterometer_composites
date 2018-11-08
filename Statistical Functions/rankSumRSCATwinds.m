function pValue = rankSumRSCATwinds(bin)

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
    
    data1u=squeeze(bin.u(iInd(i),jInd(i),1:6:end));
    countTest=squeeze(bin.count(iInd(i),jInd(i),1:6:end));
    data1u(countTest==0)=[];
    
    data1v=squeeze(bin.v(iInd(i),jInd(i),1:6:end));
    data1v(countTest==0)=[];
    
    data2u=squeeze(bin.u(iInd(i),jInd(i),2:6:end));
    countTest=squeeze(bin.count(iInd(i),jInd(i),2:6:end));
    data2u(countTest==0)=[];
    
    data2v=squeeze(bin.v(iInd(i),jInd(i),2:6:end));
    data2v(countTest==0)=[];
    
    data3u=squeeze(bin.u(iInd(i),jInd(i),3:6:end));
    countTest=squeeze(bin.count(iInd(i),jInd(i),3:6:end));
    data3u(countTest==0)=[];
    
    data3v=squeeze(bin.v(iInd(i),jInd(i),3:6:end));
    data3v(countTest==0)=[];
    
    data4u=squeeze(bin.u(iInd(i),jInd(i),4:6:end));
    countTest=squeeze(bin.count(iInd(i),jInd(i),4:6:end));
    data4u(countTest==0)=[];
    
    data4v=squeeze(bin.v(iInd(i),jInd(i),4:6:end));
    data4v(countTest==0)=[];
    
    data5u=squeeze(bin.u(iInd(i),jInd(i),5:6:end));
    countTest=squeeze(bin.count(iInd(i),jInd(i),5:6:end));
    data5u(countTest==0)=[];
    
    data5v=squeeze(bin.v(iInd(i),jInd(i),5:6:end));
    data5v(countTest==0)=[];
    
    data6u=squeeze(bin.u(iInd(i),jInd(i),6:6:end));
    countTest=squeeze(bin.count(iInd(i),jInd(i),6:6:end));
    data6u(countTest==0)=[];
    
    data6v=squeeze(bin.v(iInd(i),jInd(i),6:6:end));
    data6v(countTest==0)=[];
    
    try
        p1u=ranksum(data1u,data4u);
        p2u=ranksum(data2u,data5u);
        p3u=ranksum(data3u,data6u);
        
        p4u=ranksum(data1u,data2u);
        p5u=ranksum(data1u,data3u);
        p6u=ranksum(data1u,data5u);
        p7u=ranksum(data1u,data6u);
        
        p8u=ranksum(data2u,data3u);
        p9u=ranksum(data2u,data4u);
        p10u=ranksum(data2u,data6u);
        
        p11u=ranksum(data3u,data4u);
        p12u=ranksum(data3u,data5u);
        
        p13u=ranksum(data4u,data5u);
        p14u=ranksum(data4u,data6u);
        
        p15u=ranksum(data5u,data6u);
        
        p1v=ranksum(data1v,data4v);
        p2v=ranksum(data2v,data5v);
        p3v=ranksum(data3v,data6v);
        
        p4v=ranksum(data1v,data2v);
        p5v=ranksum(data1v,data3v);
        p6v=ranksum(data1v,data5v);
        p7v=ranksum(data1v,data6v);
        
        p8v=ranksum(data2v,data3v);
        p9v=ranksum(data2v,data4v);
        p10v=ranksum(data2v,data6v);
        
        p11v=ranksum(data3v,data4v);
        p12v=ranksum(data3v,data3v);
        
        p13v=ranksum(data4v,data5v);
        p14v=ranksum(data4v,data6v);
        
        p15v=ranksum(data5v,data6v);
        
        p=min([p1u,p2u,p3u,p4u,p5u,p6u,p7u,p8u,p9u,p10u,p11u,p12u,p13u,p14u,p15u,...
            p1v,p2v,p3v,p4v,p5v,p6v,p7v,p8v,p9v,p10v,p11v,p12v,p13v,p14v,p15v,]);
        pValue(iInd(i),jInd(i))=p;
    catch
        pValue(iInd(i),jInd(i))=nan;
        failureCount=failureCount+1;
    end
    
end

end

