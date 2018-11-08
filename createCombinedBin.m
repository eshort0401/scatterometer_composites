function bin = createCombinedBin(HY2SCATbinFileName,ASCATbinFileName,OSCATbinFileName)

binTemp=load(HY2SCATbinFileName,'bin');
binTemp=binTemp.bin;

bin.numDays=binTemp.numDays;
bin.numDaysTot=binTemp.numDaysTot;
bin.x=binTemp.x;
bin.y=binTemp.y;
bin.NtPerDay=6;
bin.binArgs=binTemp.binArgs;
bin.binArgs.dataSource='Combined 25 km';
bin.dt=0;

bin.dateCell=binTemp.dateCell;
for i=1:2
    bin.dateCell{i}([4 5 6])=[0 0 0];
end

bin.u=nan(length(bin.x),length(bin.y),bin.NtPerDay*bin.numDaysTot);
bin.v=nan(length(bin.x),length(bin.y),bin.NtPerDay*bin.numDaysTot);
bin.count=nan(length(bin.x),length(bin.y),bin.NtPerDay*bin.numDaysTot);
bin.time=nan(length(bin.x),length(bin.y),bin.NtPerDay*bin.numDaysTot);

for i=0:bin.numDaysTot-1

    bin.u(:,:,i*bin.NtPerDay+1)=binTemp.u(:,:,i*binTemp.NtPerDay+1);
    bin.u(:,:,i*bin.NtPerDay+4)=binTemp.u(:,:,i*binTemp.NtPerDay+binTemp.NtPerDay/2+1);
    
    bin.v(:,:,i*bin.NtPerDay+1)=binTemp.v(:,:,i*binTemp.NtPerDay+1);
    bin.v(:,:,i*bin.NtPerDay+4)=binTemp.v(:,:,i*binTemp.NtPerDay+binTemp.NtPerDay/2+1);
    
    bin.time(:,:,i*bin.NtPerDay+1)=binTemp.time(:,:,i*binTemp.NtPerDay+1);
    bin.time(:,:,i*bin.NtPerDay+4)=binTemp.time(:,:,i*binTemp.NtPerDay+binTemp.NtPerDay/2+1);
    
    bin.count(:,:,i*bin.NtPerDay+1)=binTemp.count(:,:,i*binTemp.NtPerDay+1);
    bin.count(:,:,i*bin.NtPerDay+4)=binTemp.count(:,:,i*binTemp.NtPerDay+binTemp.NtPerDay/2+1);

end

clear binTemp;

binTemp=load(ASCATbinFileName,'bin');
binTemp=binTemp.bin;

for i=0:bin.numDaysTot-1

    bin.u(:,:,i*bin.NtPerDay+2)=binTemp.u(:,:,i*binTemp.NtPerDay+1);
    bin.u(:,:,i*bin.NtPerDay+5)=binTemp.u(:,:,i*binTemp.NtPerDay+binTemp.NtPerDay/2+1);
    
    bin.v(:,:,i*bin.NtPerDay+2)=binTemp.v(:,:,i*binTemp.NtPerDay+1);
    bin.v(:,:,i*bin.NtPerDay+5)=binTemp.v(:,:,i*binTemp.NtPerDay+binTemp.NtPerDay/2+1);
    
    bin.time(:,:,i*bin.NtPerDay+2)=binTemp.time(:,:,i*binTemp.NtPerDay+1);
    bin.time(:,:,i*bin.NtPerDay+5)=binTemp.time(:,:,i*binTemp.NtPerDay+binTemp.NtPerDay/2+1);
    
    bin.count(:,:,i*bin.NtPerDay+2)=binTemp.count(:,:,i*binTemp.NtPerDay+1);
    bin.count(:,:,i*bin.NtPerDay+5)=binTemp.count(:,:,i*binTemp.NtPerDay+binTemp.NtPerDay/2+1);

end

clear binTemp;

binTemp=load(OSCATbinFileName,'bin');
binTemp=binTemp.bin;

for i=0:bin.numDaysTot-1

    bin.u(:,:,i*bin.NtPerDay+3)=binTemp.u(:,:,i*binTemp.NtPerDay+1);
    bin.u(:,:,i*bin.NtPerDay+6)=binTemp.u(:,:,i*binTemp.NtPerDay+binTemp.NtPerDay/2+1);
    
    bin.v(:,:,i*bin.NtPerDay+3)=binTemp.v(:,:,i*binTemp.NtPerDay+1);
    bin.v(:,:,i*bin.NtPerDay+6)=binTemp.v(:,:,i*binTemp.NtPerDay+binTemp.NtPerDay/2+1);
    
    bin.time(:,:,i*bin.NtPerDay+3)=binTemp.time(:,:,i*binTemp.NtPerDay+1);
    bin.time(:,:,i*bin.NtPerDay+6)=binTemp.time(:,:,i*binTemp.NtPerDay+binTemp.NtPerDay/2+1);
    
    bin.count(:,:,i*bin.NtPerDay+3)=binTemp.count(:,:,i*binTemp.NtPerDay+1);
    bin.count(:,:,i*bin.NtPerDay+6)=binTemp.count(:,:,i*binTemp.NtPerDay+binTemp.NtPerDay/2+1);

end


