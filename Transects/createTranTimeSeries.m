function tranTimeSeries = createTranTimeSeries(tran,times)

tranTimeSeries.distance=tran{1}.distance;
tranTimeSeries.pProj=tran{1}.pProj;
tranTimeSeries.timeSeries=zeros(length(tran{1}.distance),length(tran));
tranTimeSeries.times=times;
tranTimeSeries.label=tran{1}.label;

for i=1:length(tran)
    tranTimeSeries.timeSeries(:,i)=tran{i}.pertProj;
end

% Or simpler to use function handles?
tranTimeSeries.harmFun=cell(1,length(tran{1}.distance));

for i=1:length(tran{1}.distance)
    
    % Initialise coefficient matrix
    E=zeros(4);
    
    % Time is measured in hours here!
    EhtCos=@(h,t) cos(h*(2*pi/24)*t);
    EhtSin=@(h,t) sin(h*(2*pi/24)*t);
    
    for j=1:length(times)
        E(j,:)=[EhtCos(1,times(j)) EhtSin(1,times(j)) ...
            EhtCos(2,times(j)) EhtSin(2,times(j))];
    end
    
    % Solve E*coef=timeSeries(i,:), coef=inv(E)*timeSeries(i,:)
    
    % Note i is now counting distance!
    coef=E\tranTimeSeries.timeSeries(i,:)';
    
    tranTimeSeries.harmFun{i}=@(t) coef(1)*EhtCos(1,t)+coef(2)*EhtSin(1,t)+...
        coef(3)*EhtCos(2,t)+coef(4)*EhtSin(2,t);

end

end

