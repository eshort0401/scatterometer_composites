countRanAlt=sum(bin.count(:,:,4320-(2:17))>0,3);

uRanAv=sum(bin.u(:,:,4320-(2:17)),3);
uRanAv=uRanAv./countRanAlt;

vRanAv=sum(bin.v(:,:,4320-(2:17)),3);
vRanAv=vRanAv./countRanAlt;

for k=4320-(1:7)

    uPertK=bin.u(:,:,k)-uRanAv;
    vPertK=bin.v(:,:,k)-vRanAv;

    uPertK(countRanAlt<avArgs.runThreshAlt)=0;
    vPertK(countRanAlt<avArgs.runThreshAlt)=0;

    uPert(:,:,k)=uPertK;
    vPert(:,:,k)=vPertK;

end