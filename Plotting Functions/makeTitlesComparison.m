function [plotTitle, figFileName]=makeTitlesComparison(timeBin,startTime,endTime,dataSource1,dataSource2,timeScale,type,dt,MJOstring)

    inHours=floor((startTime(4)*60*60+startTime(5)*60+startTime(6)+(timeBin-1)*dt)/(60*60));
    inMinutes=floor((startTime(4)*60*60+startTime(5)*60+startTime(6)+(timeBin-1)*dt-inHours*60*60)/60);
    inHours=floor(rem(startTime(4)*60*60+startTime(5)*60+startTime(6)+(timeBin-1)*dt,24*60*60)/(60*60));

    fHours=floor((startTime(4)*60*60+startTime(5)*60+timeBin*dt)/(60*60));
    fMinutes=floor((startTime(4)*60*60+startTime(5)*60+startTime(6)+timeBin*dt-fHours*60*60)/60);
    fHours=floor(rem(startTime(4)*60*60+startTime(5)*60+timeBin*dt,24*60*60)/(60*60));

    plotTitle=sprintf([dataSource1 ' minus ' dataSource2 ' ' type ' between ' ...
        datestr(startTime,0) ' and ' datestr(endTime,0) ' ' timeScale ...
        '.\n %02.0f:%02.0f to %02.0f:%02.0f, MJO phases ' MJOstring],...
        inHours,inMinutes,fHours,fMinutes);

    figFileName=sprintf([type,num2str(inHours),'-',num2str(inMinutes),timeScale]);

end


