function bin = subBin(bin,startLat,endLat,startLon,endLon)
% Form a smaller data set from a larger (allready binned) data set. 

if isfield(bin,'uComp')
    
    bin.uComp=bin.uComp(bin.x>=startLon & bin.x<=endLon,bin.y>=startLat & bin.y<=endLat,:);
    bin.vComp=bin.vComp(bin.x>=startLon & bin.x<=endLon,bin.y>=startLat & bin.y<=endLat,:);
    bin.countComp=bin.countComp(bin.x>=startLon & bin.x<=endLon,bin.y>=startLat & bin.y<=endLat,:);
    bin.tComp=bin.tComp(bin.x>=startLon & bin.x<=endLon,bin.y>=startLat & bin.y<=endLat,:);

end

if isfield(bin,'uPertComp')
    
    bin.uPertComp=bin.uPertComp(bin.x>=startLon & bin.x<=endLon,bin.y>=startLat & bin.y<=endLat,:);
    bin.vPertComp=bin.vPertComp(bin.x>=startLon & bin.x<=endLon,bin.y>=startLat & bin.y<=endLat,:);

end

if isfield(bin','pValue')
    bin.pValue=bin.pValue(bin.x>=startLon & bin.x<=endLon,bin.y>=startLat & bin.y<=endLat);
end

if isfield(bin','div')
    bin.div=bin.div(bin.x>=startLon & bin.x<=endLon,bin.y>=startLat & bin.y<=endLat,:);
end

if isfield(bin,'uVar')
    bin.uVar=bin.uVar(bin.x>=startLon & bin.x<=endLon,bin.y>=startLat & bin.y<=endLat,:);
    bin.vVar=bin.vVar(bin.x>=startLon & bin.x<=endLon,bin.y>=startLat & bin.y<=endLat,:);
    bin.coVar=bin.coVar(bin.x>=startLon & bin.x<=endLon,bin.y>=startLat & bin.y<=endLat,:);
    bin.tVar=bin.tVar(bin.x>=startLon & bin.x<=endLon,bin.y>=startLat & bin.y<=endLat,:);
end

if isfield(bin,'diurnal')
    bin.diurnal.uCoef=bin.diurnal.uCoef(bin.x>=startLon & bin.x<=endLon,bin.y>=startLat & bin.y<=endLat,:);
    bin.diurnal.vCoef=bin.diurnal.vCoef(bin.x>=startLon & bin.x<=endLon,bin.y>=startLat & bin.y<=endLat,:);
    bin.diurnal.uR=bin.diurnal.uR(bin.x>=startLon & bin.x<=endLon,bin.y>=startLat & bin.y<=endLat);
    bin.diurnal.vR=bin.diurnal.vR(bin.x>=startLon & bin.x<=endLon,bin.y>=startLat & bin.y<=endLat);
    bin.diurnal.tMax=bin.diurnal.tMax(bin.x>=startLon & bin.x<=endLon,bin.y>=startLat & bin.y<=endLat);

end

if isfield(bin,'u')
    bin.u=bin.u(bin.x>=startLon & bin.x<=endLon,bin.y>=startLat & bin.y<=endLat,:);
    bin.v=bin.v(bin.x>=startLon & bin.x<=endLon,bin.y>=startLat & bin.y<=endLat,:);
    bin.time=bin.time(bin.x>=startLon & bin.x<=endLon,bin.y>=startLat & bin.y<=endLat,:);
    bin.count=bin.count(bin.x>=startLon & bin.x<=endLon,bin.y>=startLat & bin.y<=endLat,:);
end

if isfield(bin,'pComp')
    
    bin.pComp=bin.pComp(bin.x>=startLon & bin.x<=endLon,bin.y>=startLat & bin.y<=endLat,:);
%     bin.pPertComp=bin.pPertComp(bin.x>=startLon & bin.x<=endLon,bin.y>=startLat & bin.y<=endLat,:);

end

bin.x=bin.x(bin.x>=startLon & bin.x<=endLon);
bin.y=bin.y(bin.y>=startLat & bin.y<=endLat);