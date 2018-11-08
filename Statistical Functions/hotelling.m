function pValue = hotelling(bin1,bin2,time1,time2)
% Note 0 corresponds to standard hotelling, 1 corresponds to Yao
% approximation, 2 corresponds to Johansen's approx

method=2;

[iInd, jInd] = find(bin1.uPertVar(:,:,time1)~=0 &...
        bin2.vPertVar(:,:,time2)~=0);
    
T2=nan(length(bin1.x),length(bin1.y));
nu=nan(length(bin1.x),length(bin1.y));

for m=1:length(iInd)

    S1=zeros(2,2);
    S1(1,1)=bin1.uPertVar(iInd(m),jInd(m),time1);
    S1(2,2)=bin1.vPertVar(iInd(m),jInd(m),time1);
    S1(1,2)=bin1.pertCoVar(iInd(m),jInd(m),time1);
    S1(2,1)=S1(1,2);

    S2=zeros(2,2);
    S2(1,1)=bin2.uPertVar(iInd(m),jInd(m),time2);
    S2(2,2)=bin2.vPertVar(iInd(m),jInd(m),time2);
    S2(1,2)=bin2.pertCoVar(iInd(m),jInd(m),time2);
    S2(2,1)=S2(1,2);

    n1=bin1.countComp(iInd(m),jInd(m),time1);
    n2=bin2.countComp(iInd(m),jInd(m),time2);

    y1=zeros(2,1);
    y1(1)=bin1.uPertComp(iInd(m),jInd(m),time1);
    y1(2)=bin1.vPertComp(iInd(m),jInd(m),time1);

    y2=zeros(2,1);
    y2(1)=bin2.uPertComp(iInd(m),jInd(m),time2);
    y2(2)=bin2.vPertComp(iInd(m),jInd(m),time2);
    
    if method==0 % Standard Hotelling
        Sp=((n1-1).*S1+(n2-1).*S2)/(n1+n2-2);
        T2(iInd(m),jInd(m))=((n1*n2)/(n1+n2))*(y1-y2)'*(Sp\(y1-y2));
    elseif method==1 % Yao
        V1=S1./n1;
        V2=S2./n2;
        Se=V1+V2;
        % Calculate T*2
        T2(iInd(m),jInd(m))=(y1-y2)'*(Se\(y1-y2));

        nuInv1=(1/(n1-1)).*(((y1-y2)'*(Se\V1)*(Se\(y1-y2)))^2);
        nuInv2=(1/(n2-1)).*(((y1-y2)'*(Se\V2)*(Se\(y1-y2)))^2);

        nuInv=(1/T2(iInd(m),jInd(m)))*(nuInv1+nuInv2);
        nu(iInd(m),jInd(m))=1/nuInv;
    elseif method==2 % Merwe
        V1=S1./n1;
        V2=S2./n2;
        Se=V1+V2;
        
        nuNum=Se(1,1)^2+Se(2,2)^2+trace(Se)^2;
        nuDen=(1/(n1-1))*(V1(1,1)^2+V1(2,2)^2+trace(V1)^2)+...
            (1/(n2-1))*(V2(1,1)^2+V2(2,2)^2+trace(V2)^2);
        nu(iInd(m),jInd(m))=nuNum/nuDen;
        T2(iInd(m),jInd(m))=(y1-y2)'*(Se\(y1-y2));
    end
end

if method==0 % Standard Hotelling
    F=T2.*(bin1.countComp(:,:,time1)+bin2.countComp(:,:,time2)-2-1)./...
        ((bin1.countComp(:,:,time1)+bin2.countComp(:,:,time2)-2)*2);
    nu1=2;
    nu2=bin1.countComp(:,:,time1)+bin2.countComp(:,:,time2)-2-1;
    pValue=1-fcdf(F,nu1,nu2);
elseif method==1 % Yao's method
    p=2;
    F=((nu-p+1)./(nu*p)).*T2;
    pValue=1-fcdf(F,p,nu-p+1); 
elseif method==2 % Merwe method
    p=2;
    F=((nu-p+1)./(nu*p)).*T2;
    pValue=1-fcdf(F,p,nu-p+1); 
end


end

