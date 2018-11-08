function pValue = rankSumERA(tp,lat,lon)

failureCount=0;

pValue=nan(length(lon),length(lat));

for i=1:size(tp,1)
    for j=1:size(tp,2)
    
    % Now were only working with the spatial point (iInd(i),jInd(i)).
    % First create array of the form (u,v,t) for this point.
    
        data1=squeeze(tp(i,j,1:8:end));
        data2=squeeze(tp(i,j,3:8:end));
        data3=squeeze(tp(i,j,5:8:end));
        data4=squeeze(tp(i,j,7:8:end));

        try
            p1=ranksum(data1,data2);
            p2=ranksum(data1,data3);
            p3=ranksum(data1,data4);
            
            p4=ranksum(data2,data3);
            p5=ranksum(data2,data4);
            
            p6=ranksum(data3,data4);

            p=min([p1,p2,p3,p4,p5,p6]);
            pValue(i,j)=p;
        catch
            pValue(i,j)=nan;
            failureCount=failureCount+1;
        end
    end
end

end

