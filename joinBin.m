function bin = joinBin(bin1,bin2)
% Combine two data sets. Note bin1.x must continue on to bin2.x, and the
% two y arrays must be equal.

if ~isfield(bin1,'dt')
    bin1.dt=24*3600/bin1.NtPerDay;
    bin2.dt=24*3600/bin1.NtPerDay;
    bin.dt=24*3600/bin1.NtPerDay;
end

dLon=round(bin1.x(2)-bin1.x(1),8);
if (bin2.x(1)~=bin1.x(end)+dLon || any(bin1.y ~=bin2.y) || ...
        bin1.dt ~= bin2.dt)
    error('Binned data sets incompatible.');
end

bin.dateCell=bin1.dateCell;
bin.numDays=bin1.numDays;
bin.numDaysTot=bin1.numDaysTot;
bin.dt=bin1.dt;
bin.x=[bin1.x,bin2.x];
bin.y=bin1.y;
bin.binArgs=bin1.binArgs;
bin.avArgs=bin1.avArgs;
bin.NtPerDay=bin1.NtPerDay;

fprintf('Joining bins. \n');

if isfield(bin1,'uComp')
    
    bin.uComp=cat(1,bin1.uComp,bin2.uComp);
    bin.vComp=cat(1,bin1.vComp,bin2.vComp);
    bin.countComp=cat(1,bin1.countComp,bin2.countComp);
    bin.tComp=cat(1,bin1.tComp,bin2.tComp);

end

if isfield(bin1,'pComp')
    bin.pComp=cat(1,bin1.pComp,bin2.pComp);
    bin.countComp=cat(1,bin1.countComp,bin2.countComp);
    bin.tComp=cat(1,bin1.tComp,bin2.tComp);
    bin.pVar=cat(1,bin1.pVar,bin2.pVar);
end

if isfield(bin1,'pPertComp')
    bin.pPertComp=cat(1,bin1.pPertComp,bin2.pPertComp);
    bin.pPertVar=cat(1,bin1.pPertVar,bin2.pPertVar);
end
    
if isfield(bin1,'uPertComp')
    
    bin.uPertComp=cat(1,bin1.uPertComp,bin2.uPertComp);
    bin.vPertComp=cat(1,bin1.vPertComp,bin2.vPertComp);

end

if isfield(bin1,'pValue')
    bin.pValue=cat(1,bin1.pValue,bin2.pValue);
end

if isfield(bin1,'uVar')
    bin.uVar=cat(1,bin1.uVar,bin2.uVar);
    bin.vVar=cat(1,bin1.vVar,bin2.vVar);
    bin.coVar=cat(1,bin1.coVar,bin2.coVar);
    bin.tVar=cat(1,bin1.tVar,bin2.tVar);
end

if isfield(bin1,'uPertVar')
    bin.uPertVar=cat(1,bin1.uPertVar,bin2.uPertVar);
    bin.vPertVar=cat(1,bin1.vPertVar,bin2.vPertVar);
    bin.pertCoVar=cat(1,bin1.pertCoVar,bin2.pertCoVar);
end

if isfield(bin1,'pVar')
    bin.pVar=cat(1,bin1.pVar,bin2.pVar);
end

end

