function [plotTitle, figFileName]=makeTitles(timeBin,startTime,endTime,dataSource,timeScale,type,dt)

    inHours=floor((startTime(4)*60*60+startTime(5)*60+startTime(6)+(timeBin-1)*dt)/(60*60));
    inMinutes=floor((startTime(4)*60*60+startTime(5)*60+startTime(6)+(timeBin-1)*dt-inHours*60*60)/60);
    inHours=floor(rem(startTime(4)*60*60+startTime(5)*60+startTime(6)+(timeBin-1)*dt,24*60*60)/(60*60));

    fHours=floor((startTime(4)*60*60+startTime(5)*60+timeBin*dt)/(60*60));
    fMinutes=floor((startTime(4)*60*60+startTime(5)*60+startTime(6)+timeBin*dt-fHours*60*60)/60);
    fHours=floor(rem(startTime(4)*60*60+startTime(5)*60+timeBin*dt,24*60*60)/(60*60));
    
    plotTitle=sprintf('%02d%02d-%02d%02d',...
    inHours,inMinutes,fHours,fMinutes);
    
    figFileName=sprintf([type,plotTitle]);

end


