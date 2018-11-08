function bin = diurnalCycle(binCell)
% Solve for diurnal harmonics

kInd=cell(length(binCell),1);

% Assume each bin has only two nonzero time components
diurnal.uCoef=nan(length(binCell{1}.x),length(binCell{1}.y),3);
diurnal.vCoef=nan(length(binCell{1}.x),length(binCell{1}.y),3);
diurnal.uR=nan(length(binCell{1}.x),length(binCell{1}.y));
diurnal.vR=nan(length(binCell{1}.x),length(binCell{1}.y));
diurnal.tMax=nan(length(binCell{1}.x),length(binCell{1}.y));

dataExists=true(length(binCell{1}.x),length(binCell{1}.y));
    
for i=1:length(binCell)
    dataExists=dataExists & any(binCell{i}.countComp>0,3); 
end

[iInd, jInd] = find(dataExists);

% Iterate through grid cells
for m=1:length(iInd)
    
    % Find nonzero time bins.
    numObs=zeros(1,length(binCell));
    cumObs=zeros(1,length(binCell));
    totalObs=0;
    for i=1:length(binCell)
        kInd{i}=find(squeeze(squeeze(binCell{i}.countComp(iInd(m),jInd(m),:)>0)));
        kInd{i}=unique(kInd{i});
        numObs(i)=length(kInd{i});
        cumObs(i)=totalObs;
        totalObs=totalObs+numObs(i);
    end
    
    % Initialise vectors for regression
    E=ones(totalObs,3);
    u=nan(totalObs,1);
    v=nan(totalObs,1);
    t=nan(totalObs,1);

    % Iterate through data sets
    for i=1:length(binCell)
        
        % Iterate through time cells of data set i.
        for j=1:length(kInd{i})
        
            obsNum=cumObs(i)+j;
            
            E(obsNum,2)=cos((2*pi/24)*binCell{i}.tComp(iInd(m),jInd(m),kInd{i}(j)))';
            E(obsNum,3)=sin((2*pi/24)*binCell{i}.tComp(iInd(m),jInd(m),kInd{i}(j)))';

            u(obsNum)=binCell{i}.uComp(iInd(m),jInd(m),kInd{i}(j));
            v(obsNum)=binCell{i}.vComp(iInd(m),jInd(m),kInd{i}(j));
            t(obsNum)=binCell{i}.tComp(iInd(m),jInd(m),kInd{i}(j));
            
        end
    
    end
    
    if (any(isnan(E(:))) || any(isnan(u)) || any(isnan(v)) || any(isnan(t)))
        continue;
    end
    
    diurnal.uCoef(iInd(m),jInd(m),:)=E\u;
    diurnal.vCoef(iInd(m),jInd(m),:)=E\v;
    
    % Calculate R^2 coefficient
    uSStot=sum((u-mean(u)).^2);
    vSStot=sum((v-mean(v)).^2);
    
    uModelled=diurnal.uCoef(iInd(m),jInd(m),1)+...
        diurnal.uCoef(iInd(m),jInd(m),2)*cos((2*pi/24)*t')+...
        diurnal.uCoef(iInd(m),jInd(m),3)*sin((2*pi/24)*t');
    
    vModelled=diurnal.vCoef(iInd(m),jInd(m),1)+...
        diurnal.vCoef(iInd(m),jInd(m),2)*cos((2*pi/24)*t')+...
        diurnal.vCoef(iInd(m),jInd(m),3)*sin((2*pi/24)*t');
    
    uSSres=sum((uModelled'-u).^2);
    vSSres=sum((vModelled'-v).^2);
    
    diurnal.uR(iInd(m),jInd(m))=1-uSSres/uSStot;
    diurnal.vR(iInd(m),jInd(m))=1-vSSres/vSStot;
    
    diurnal.tMax(iInd(m),jInd(m))=mod((6/pi)*atan(2*(...
        diurnal.uCoef(iInd(m),jInd(m),2)*diurnal.uCoef(iInd(m),jInd(m),3)+...
        diurnal.vCoef(iInd(m),jInd(m),2)*diurnal.vCoef(iInd(m),jInd(m),3))/(...
        diurnal.uCoef(iInd(m),jInd(m),2)^2+diurnal.vCoef(iInd(m),jInd(m),2)^2-...
        diurnal.uCoef(iInd(m),jInd(m),3)^2-diurnal.vCoef(iInd(m),jInd(m),3)^2)),12);
    
    secondDeriv=-diurnal.uCoef(iInd(m),jInd(m),2)*diurnal.uCoef(iInd(m),jInd(m),3)*sin(2*pi*diurnal.tMax(iInd(m),jInd(m))/12)-...
    diurnal.vCoef(iInd(m),jInd(m),2)*diurnal.vCoef(iInd(m),jInd(m),3)*sin(2*pi*diurnal.tMax(iInd(m),jInd(m))/12);

    if secondDeriv>0 
        diurnal.tMax(iInd(m),jInd(m))=mod(diurnal.tMax(iInd(m),jInd(m))+6,12);
    elseif secondDeriv==0
        diurnal.tMax(iInd(m),jInd(m))=NaN;
    end
end

bin=binCell{1};
bin.diurnal=diurnal;

end

