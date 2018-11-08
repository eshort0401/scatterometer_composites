function [distance, compProj, pertProj, pertVarProj, pProj] = ...
    tranProj(binAv,startLat,startLon,endLat,endLon,nPoints,time)
% tranProj ----------------------------------------------------------
%
% INPUT: points is the number of points in the transect.
%
% OUTPUT: proj is a 1D array containing the projection.
%--------------------------------------------------------------------------

dx=(endLon-startLon)/(nPoints-1);
dy=(endLat-startLat)/(nPoints-1);

[X, Y]=meshgrid(binAv.x,binAv.y);
[interpX, interpY]=meshgrid(startLon:dx:endLon, startLat:dy:endLat);

% Use approximations from wiki.
distance=(0:sqrt((dx*110.574)^2+(dy*111.320)^2):sqrt((dx*110.574)^2+(dy*111.320)^2)*(nPoints-1));

if isfield(binAv,'uPertComp')

    binAv.uPertComp(binAv.uPertComp==0)=NaN;
    binAv.vPertComp(binAv.vPertComp==0)=NaN;
    binAv.uComp(binAv.uComp==0)=NaN;
    binAv.vComp(binAv.vComp==0)=NaN;
    binAv.uVar(binAv.uVar==0)=NaN;
    binAv.vVar(binAv.vVar==0)=NaN;
    binAv.coVar(binAv.coVar==0)=NaN;
    
    uCompProj = diag(interp2(X,Y,binAv.uComp(:,:,time)',interpX,interpY));
    vCompProj = diag(interp2(X,Y,binAv.vComp(:,:,time)',interpX,interpY));
    uPertProj = diag(interp2(X,Y,binAv.uPertComp(:,:,time)',interpX,interpY));
    vPertProj = diag(interp2(X,Y,binAv.vPertComp(:,:,time)',interpX,interpY));
    uPertVarProj = diag(interp2(X,Y,binAv.uPertVar(:,:,time)',interpX,interpY));
    pertCoVarProj = diag(interp2(X,Y,binAv.pertCoVar(:,:,time)',interpX,interpY));
    pProj = diag(interp2(X,Y,binAv.pValue(:,:)',interpX,interpY));

    compProj=zeros(1,nPoints);
    pertProj=zeros(1,nPoints);
    pertVarProj=zeros(1,nPoints);

    for i=1:nPoints
        compProj(i)=([uCompProj(i) vCompProj(i)]*[dx;dy])./sqrt(dx^2+dy^2);
        pertProj(i)=([uPertProj(i) vPertProj(i)]*[dx;dy])./sqrt(dx^2+dy^2);
        pertVarProj(i)=(1/sqrt(dx^2+dy^2)).*(dx*uPertVarProj(i)+dy*pertCoVarProj(i));
    end
    
elseif isfield(binAv,'pPertComp')
    
    binAv.pPertComp(binAv.pPertComp==0)=NaN;
    binAv.pComp(binAv.pComp==0)=NaN;
    binAv.pVar(binAv.pVar==0)=NaN;
    
    compProj = diag(interp2(X,Y,binAv.pComp(:,:,time)',interpX,interpY));
    pertProj = diag(interp2(X,Y,binAv.pPertComp(:,:,time)',interpX,interpY));
    pertVarProj = diag(interp2(X,Y,binAv.pPertVar(:,:,time)',interpX,interpY));
    pProj = diag(interp2(X,Y,binAv.pValue(:,:)',interpX,interpY));
    
end





