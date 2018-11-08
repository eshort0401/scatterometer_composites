function [statSig, pStar] = FDRsignificance(pValue,alphaFDR)
% Determine where data is statistically significant using the specified
% alpha, alphaFDR and the FDR procedure of Wilks. statSig is an array of
% same dimension as pValue with ones if data is statistically significant,
% 0 otherwise. 

pList=pValue(:);
pList(isnan(pList) | isinf(pList))=[];
pList=sort(pList);

N=length(pList);

pStar=alphaFDR.*(1:N)'./N;
pStar=find(pList<=pStar,1,'last');
pStar=pList(pStar);

statSig=pValue<=pStar;
end

