function bin = concatenateBin(bin1,bin2)
% Combine two data sets. Note bin1.x must continue on to bin2.x, and the
% two y arrays must be equal.

if ~isfield(bin1,'dt')
    bin1.dt=24*3600/bin1.NtPerDay;
    bin2.dt=24*3600/bin1.NtPerDay;
    bin.dt=24*3600/bin1.NtPerDay;
end

if (any(bin2.x~=bin2.x) || any(bin1.y~=bin2.y) || ...
        bin1.dt ~= bin2.dt)
    error('Binned data sets incompatible.');
end

bin.dateCell=[bin1.dateCell bin2.dateCell];
bin.numDays=[bin1.numDays bin2.numDays];
bin.numDaysTot=bin1.numDaysTot+bin2.numDaysTot;
bin.dt=bin1.dt;
bin.x=bin1.x;
bin.y=bin1.y;
bin.binArgs=bin1.binArgs;
bin.NtPerDay=bin1.NtPerDay;

fprintf('Joining bins. \n');

if isfield(bin1,'p')
    bin.p=cat(3,bin1.p,bin2.p);
    bin.count=cat(3,bin1.count,bin2.count);
    bin.time=cat(3,bin1.time,bin2.time);
end

end

