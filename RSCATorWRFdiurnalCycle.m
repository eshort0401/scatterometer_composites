function bin = RSCATorWRFdiurnalCycle(bin)

% Solve for the RSCAT diurnal cycle
E=ones(bin.NtPerDay,3);
diurnal.uCoef=nan(length(bin.x),length(bin.y),3);
diurnal.vCoef=nan(length(bin.x),length(bin.y),3);
diurnal.uR=nan(length(bin.x),length(bin.y));
diurnal.vR=nan(length(bin.x),length(bin.y));
diurnal.tMax=nan(length(bin.x),length(bin.y));

% Find grid cells with positive count
[iInd, jInd] = find(~isnan(bin.uComp(:,:,1)) &...
        ~isnan(bin.vComp(:,:,1)) & ~isnan(bin.tComp(:,:,1)));

% For each grid cell
for m=1:length(iInd)

    % Form equations and solve least squares
    E(:,2)=cos((2*pi/24)*squeeze(bin.tComp(iInd(m),jInd(m),:)))';
    E(:,3)=sin((2*pi/24)*squeeze(bin.tComp(iInd(m),jInd(m),:)))';
    
    u=squeeze(squeeze(bin.uComp(iInd(m),jInd(m),:)));
    v=squeeze(squeeze(bin.vComp(iInd(m),jInd(m),:)));
    
    if (any(isnan(E(:))) || any(isnan(u)) || any(isnan(v)))
        continue;
    end
    
    diurnal.uCoef(iInd(m),jInd(m),:)=E\u;
    diurnal.vCoef(iInd(m),jInd(m),:)=E\v;
    
    % Calculate R^2 coefficient
    uSStot=sum((u-mean(u)).^2);
    vSStot=sum((v-mean(v)).^2);
    
    uModelled=diurnal.uCoef(iInd(m),jInd(m),1)+...
        diurnal.uCoef(iInd(m),jInd(m),2)*cos((2*pi/24)*squeeze(bin.tComp(iInd(m),jInd(m),:)))'+...
        diurnal.uCoef(iInd(m),jInd(m),3)*sin((2*pi/24)*squeeze(bin.tComp(iInd(m),jInd(m),:)))';
    
    vModelled=diurnal.vCoef(iInd(m),jInd(m),1)+...
        diurnal.vCoef(iInd(m),jInd(m),2)*cos((2*pi/24)*squeeze(bin.tComp(iInd(m),jInd(m),:)))'+...
        diurnal.vCoef(iInd(m),jInd(m),3)*sin((2*pi/24)*squeeze(bin.tComp(iInd(m),jInd(m),:)))';
    
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

bin.diurnal=diurnal;

end

