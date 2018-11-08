function tran = avTrans(...
    binAv,startLat,startLon,endLatCoast,endLonCoast,...
    endLonTran,time,label)
% avTran ----------------------------------------------------------
%
% Generate multiple transects perpendicular to coast line and take averages
% to form mean transect. 
%
% OUTPUT: Average of multiple transects.
%--------------------------------------------------------------------------

% First solve linear system to determine equation for coast
A=[startLon 1;endLonCoast 1];
B=[startLat; endLatCoast];
X=A\B;
m=X(1);
cCoast=X(2);
clear A B X;
cTran=startLat+(1/m)*startLon;
endLatTran=(-1/m)*endLonTran+cTran;
if endLatTran>6 || endLatTran<-13
    error('endLatTran out of domain, change endLonTran')
end
% Calculate nTrans and nPoints
dCoast=sqrt((startLat-endLatCoast)^2+(startLon-endLonCoast)^2);
dTrans=sqrt((startLat-endLatTran)^2+(startLon-endLonTran)^2);
% Choose smaller of x and y resolutions
nTrans=floor(dCoast/min(binAv.x(2)-binAv.x(1),binAv.y(2)-binAv.y(1)));
nPoints=floor(dTrans/min(binAv.x(2)-binAv.x(1),binAv.y(2)-binAv.y(1)));

% Store compProj and pertProj vectors in cell array.
% transects=cell(1,nTrans);
iCompProj=nan(nTrans,nPoints);
iPertProj=nan(nTrans,nPoints);
iVarProj=nan(nTrans,nPoints);
iPproj=nan(nTrans,nPoints);

dxC=(endLonCoast-startLon)/(nTrans-1);
dyC=(endLatCoast-startLat)/(nTrans-1);

% dxT=(endLonTran-startLon)/(nPoints-1);
% dyT=(endLatTran-startLat)/(nPoints-1);

% [X, Y]=meshgrid(binAv.x,binAv.y);
% 
% xCoast=startLon:dxC:endLonCoast;
% yCoast=startLat:dyC:endLonCoast;

[distance, iCompProj(1,:), iPertProj(1,:), iVarProj(1,:),iPproj(1,:)] = ...
    tranProj(binAv,startLat,startLon,endLatTran,endLonTran,nPoints,time);

for i=2:nTrans
    startLonAlt=startLon+(i-1)*dxC;
    startLatAlt=startLat+(i-1)*dyC;
    endLonTranAlt=endLonTran+(i-1)*dxC;
    endLatTranAlt=endLatTran+(i-1)*dyC;
    
    [~, iCompProj(i,:), iPertProj(i,:), iVarProj(i,:),iPproj(i,:)] = ...
    tranProj(binAv,startLatAlt,startLonAlt,endLatTranAlt,endLonTranAlt,nPoints,time);
end
    
compProj=nanmean(iCompProj,1);
pertProj=nanmean(iPertProj,1);
varProj=nanmean(iVarProj,1);
pProj=nanmean(iPproj,1);

tran.distance=distance;
tran.compProj=compProj;
tran.pertProj=pertProj;
tran.varProj=varProj;
tran.pProj=pProj;

% Create data structure
tran.startLat=startLat;
tran.startLon=startLon;
tran.endLatCoast=endLatCoast;
tran.endLonCoast=endLonCoast;
tran.endLatTran=endLatTran;
tran.endLonTran=endLonTran;
tran.nTrans=nTrans;
tran.nPoints=nPoints;
tran.time=time;
tran.label=label;


end

