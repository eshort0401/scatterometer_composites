function bin = restrictRSCAT(bin)

u=zeros(length(bin.x),length(bin.y),bin.NtPerDay*bin.numDaysTot);
v=zeros(length(bin.x),length(bin.y),bin.NtPerDay*bin.numDaysTot);
count=zeros(length(bin.x),length(bin.y),bin.NtPerDay*bin.numDaysTot);
time=zeros(length(bin.x),length(bin.y),bin.NtPerDay*bin.numDaysTot);

u(:,:,1:bin.NtPerDay:end)=bin.u(:,:,1:bin.NtPerDay:end);
u(:,:,(bin.NtPerDay/2)+1:bin.NtPerDay:end)=...
    bin.u(:,:,(bin.NtPerDay/2)+1:bin.NtPerDay:end);

v(:,:,1:bin.NtPerDay:end)=bin.v(:,:,1:bin.NtPerDay:end);
v(:,:,(bin.NtPerDay/2)+1:bin.NtPerDay:end)=...
    bin.v(:,:,(bin.NtPerDay/2)+1:bin.NtPerDay:end);

count(:,:,1:bin.NtPerDay:end)=bin.count(:,:,1:bin.NtPerDay:end);
count(:,:,(bin.NtPerDay/2)+1:bin.NtPerDay:end)=...
    bin.count(:,:,(bin.NtPerDay/2)+1:bin.NtPerDay:end);

time(:,:,1:bin.NtPerDay:end)=bin.time(:,:,1:bin.NtPerDay:end);
time(:,:,(bin.NtPerDay/2)+1:bin.NtPerDay:end)=...
    bin.time(:,:,(bin.NtPerDay/2)+1:bin.NtPerDay:end);

bin.u=u;
bin.v=v;
bin.count=count;
bin.time=time;

end

