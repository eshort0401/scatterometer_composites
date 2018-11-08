function [tranJavaNorth, tranJavaSouth, ...
    tranNewGuineaNorth, tranSumatraSouthWest, ...
    tranGulfCarp, tranGulfPap, tranFlamingoBay, ...
    tranWestDarwin] = createTransects(binAv)

ind=find(squeeze(any(any(binAv.countComp,1)>0)));

tranSumatraSouthWest=cell(length(ind),1);
tranJavaSouth=cell(length(ind),1);
tranJavaNorth=cell(length(ind),1);
tranNewGuineaNorth=cell(length(ind),1);
tranGulfCarp=cell(length(ind),1);
tranGulfPap=cell(length(ind),1);
tranFlamingoBay=cell(length(ind),1);
tranWestDarwin=cell(length(ind),1);

for i=1:length(ind)
    tranSumatraSouthWest{i}=avTrans(binAv,-2,100,-6.5,105,94,ind(i),[binAv.binArgs.dataSource ' Sumatra South West' num2str(i)]);
    tranJavaSouth{i}=avTrans(binAv,-7.4,105,-8.9,114.5,104.4,ind(i),[binAv.binArgs.dataSource ' Java South' num2str(i)]);
    tranJavaNorth{i}=avTrans(binAv,-6,109,-7,114,109.6,ind(i),[binAv.binArgs.dataSource ' Java North' num2str(i)]);
    tranNewGuineaNorth{i}=avTrans(binAv,-1,137.5,-4.5,146.5,139.75,ind(i),[binAv.binArgs.dataSource ' North New Guinea ' num2str(i)]);
    tranGulfCarp{i}=avTrans(binAv,-10.51,141.88,-12.74,141.24,137.37,ind(i),[binAv.binArgs.dataSource ' Gulf Carp ' num2str(i)]);
    tranGulfPap{i}=avTrans(binAv,-8.53,144,-8.26,146,144.465,ind(i),[binAv.binArgs.dataSource ' Gulf Papua ' num2str(i)]);
    tranFlamingoBay{i}=avTrans(binAv,-5.8,138.1,-7.4,138.2,135,ind(i),[binAv.binArgs.dataSource ' Yos Sudarno ' num2str(i)]);
    tranWestDarwin{i}=avTrans(binAv,-12.8,129.93,-11.45,129.98,126,ind(i),[binAv.binArgs.dataSource ' West Darwin ' num2str(i)]);
end

end