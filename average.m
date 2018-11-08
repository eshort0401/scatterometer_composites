function bin = average(bin,avArgs)

%--------------------------------------------------------------------------
% Calculate composites and statistics
%--------------------------------------------------------------------------

if nargin<2
    avArgs.composite=1;
    avArgs.pertComp=1;
    avArgs.divergence=0;
    avArgs.t2Test=0;
    
    % Remember the perturbation variances are the same as the normal
    % variances, unless running means are used!
    avArgs.manTest=1;
    avArgs.runMean=1;
    % How many days either side to use for running mean.
    avArgs.avRan=1;
    
    % compThresh stipulates the minimum ratio a count can have to the
    % maximum count to be included in the results.
    avArgs.compThresh=0;
    % Specify minimum ratio of data points for running mean to be counted
    avArgs.runThresh=0;
    % Specify minimum number of distinct time grid cells for running mean
    % to be counted
    avArgs.runThreshAlt=0;
end

Nt=bin.numDaysTot*bin.NtPerDay;

fprintf('Calculating composites. \n');

% Calculate composites.
bin.uComp=zeros(length(bin.x),length(bin.y),bin.NtPerDay);
bin.vComp=zeros(length(bin.x),length(bin.y),bin.NtPerDay);
bin.tComp=zeros(length(bin.x),length(bin.y),bin.NtPerDay);

% bin.tVar=zeros(length(bin.x),length(bin.y),bin.NtPerDay);
bin.uVar=zeros(length(bin.x),length(bin.y),bin.NtPerDay);
bin.vVar=zeros(length(bin.x),length(bin.y),bin.NtPerDay);
bin.coVar=zeros(length(bin.x),length(bin.y),bin.NtPerDay);
bin.tVar=zeros(length(bin.x),length(bin.y),bin.NtPerDay);

% Create array to store composite count.    
bin.countComp=zeros(length(bin.x),length(bin.y),bin.NtPerDay);
if ~isfield(bin,'time')
    bin.time=zeros(length(bin.x),length(bin.y),Nt);
end

for i=1:bin.NtPerDay
    bin.uComp(:,:,i)=sum(bin.u(:,:,i+(0:bin.numDaysTot-1)*bin.NtPerDay),3);
    bin.vComp(:,:,i)=sum(bin.v(:,:,i+(0:bin.numDaysTot-1)*bin.NtPerDay),3);
    bin.tComp(:,:,i)=sum(bin.time(:,:,i+(0:bin.numDaysTot-1)*bin.NtPerDay),3);
    bin.countComp(:,:,i)=sum(bin.count(:,:,i+(0:bin.numDaysTot-1)*bin.NtPerDay)>0,3);
end

bin.uComp=bin.uComp./bin.countComp;
bin.vComp=bin.vComp./bin.countComp;
bin.tComp=bin.tComp./bin.countComp;

bin.uComp(isnan(bin.uComp) | isinf(bin.uComp))=0;
bin.vComp(isnan(bin.vComp) | isinf(bin.vComp))=0;
bin.tComp(isnan(bin.tComp) | isinf(bin.tComp))=0;

fprintf('Calculating variances. \n');

uRes=bin.u-repmat(bin.uComp,[1 1 bin.numDaysTot]);
vRes=bin.v-repmat(bin.vComp,[1 1 bin.numDaysTot]);
tRes=bin.time-repmat(bin.tComp,[1 1 bin.numDaysTot]);

% Remove residuals where we are simply subtracting binAv.uComp from zero
uRes(bin.count==0)=0;
vRes(bin.count==0)=0;
tRes(bin.count==0)=0;

% Calculate sample variance and covariance
for i=1:bin.NtPerDay
    bin.uVar(:,:,i)=sum(uRes(:,:,i+(0:bin.numDaysTot-1)*bin.NtPerDay).^2,3);
    bin.vVar(:,:,i)=sum(vRes(:,:,i+(0:bin.numDaysTot-1)*bin.NtPerDay).^2,3);
    bin.coVar(:,:,i)=sum(uRes(:,:,i+(0:bin.numDaysTot-1)*bin.NtPerDay).*...
        vRes(:,:,i+(0:bin.numDaysTot-1)*bin.NtPerDay),3);
    bin.tVar(:,:,i)=sum(tRes(:,:,i+(0:bin.numDaysTot-1)*bin.NtPerDay).^2,3);
end

bin.uVar=bin.uVar./(bin.countComp-1);
bin.vVar=bin.vVar./(bin.countComp-1);
bin.coVar=bin.coVar./(bin.countComp-1);
bin.tVar=bin.tVar./(bin.countComp-1); 

bin.uVar(isnan(bin.uVar) | isinf(bin.uVar) | bin.countComp<=1)=0;
bin.vVar(isnan(bin.vVar) | isinf(bin.vVar) | bin.countComp<=1)=0;
bin.coVar(isnan(bin.coVar) | isinf(bin.coVar) | bin.countComp<=1)=0;
bin.tVar(isnan(bin.tVar) | isinf(bin.tVar) | bin.countComp<=1)=0;

clear uRes vRes tRes

bin=rmfield(bin,'time');
startTime=bin.dateCell{1};
bin.tComp(bin.tComp==0)=NaN;
bin.tComp=mod(bin.tComp+startTime(4)+startTime(5)/60+startTime(6)/3600,24);

%--------------------------------------------------------------------------
% Calculate perturbations if required
%--------------------------------------------------------------------------
    
kInd=find(squeeze(squeeze(any(any(bin.countComp>0,1),2)))==1);
kInd=unique(kInd);
avArgs.testTime1=kInd(1);
avArgs.testTime2=kInd(2);

% Calculate perturbation winds if required. Be careful to not double weight
% the current time. For now ignore future most time. 
if avArgs.pertComp
    
    fprintf('Calculating perturbation winds. \n');
   
    bin.uPert=zeros(length(bin.x),length(bin.y),Nt);
    bin.vPert=zeros(length(bin.x),length(bin.y),Nt);
    
    % Use basic mean.  
    if avArgs.runMean==0
        
        uCompNan=bin.uComp;
        vCompNan=bin.vComp;
        
        % Check there is data at each time.
        uCompNan(bin.uComp==0)=NaN;
        vCompNan(bin.vComp==0)=NaN;

        uRanAv=nanmean(uCompNan,3);
        vRanAv=nanmean(vCompNan,3);

        clear uCompNan vCompNan
        
        for i=1:size(bin.u,3)
            bin.uPert(:,:,i)=bin.u(:,:,i)-uRanAv;
            bin.vPert(:,:,i)=bin.v(:,:,i)-vRanAv;
        end
        
        bin.uPert(isnan(bin.uPert) | isinf(bin.uPert) | bin.count==0)=0;
        bin.vPert(isnan(bin.uPert) | isinf(bin.uPert) | bin.count==0)=0;
       
    % Use simple running mean method.
    elseif avArgs.runMean==1
        % Check running mean not too long.
        if avArgs.avRan*2>=Nt
            error('Running mean longer than number of time cells');
        end
        
        bin.uPert=zeros(length(bin.x),length(bin.y),Nt);
        bin.vPert=zeros(length(bin.x),length(bin.y),Nt);

        for k=avArgs.avRan*bin.NtPerDay+1:Nt-avArgs.avRan*bin.NtPerDay+1
            
            countRanAlt=sum(bin.count(:,:,k-avArgs.avRan*bin.NtPerDay:k+avArgs.avRan*bin.NtPerDay-1)>0,3);

            uRanAv=sum(bin.u(:,:,k-avArgs.avRan*bin.NtPerDay:k+avArgs.avRan*bin.NtPerDay-1),3);
            uRanAv=uRanAv./countRanAlt;
            
            vRanAv=sum(bin.v(:,:,k-avArgs.avRan*bin.NtPerDay:k+avArgs.avRan*bin.NtPerDay-1),3);
            vRanAv=vRanAv./countRanAlt;
            
            bin.uPert(:,:,k)=bin.u(:,:,k)-uRanAv;
            bin.vPert(:,:,k)=bin.v(:,:,k)-vRanAv;
            
        end
        
        bin.uPert(isinf(bin.uPert) | isnan(bin.uPert) | bin.count==0)=0;
        bin.vPert(isinf(bin.uPert) | isnan(bin.uPert) | bin.count==0)=0;
        
%         % Remove runMean values around the ends of disjoint intervals 
%         for i=1:length(bin.numDays)-1
%             intervalBoundary=sum(bin.numDays(1:i))*bin.NtPerDay-...
%                 avArgs.avRan*bin.NtPerDay:sum(bin.numDays(1:i))*...
%                 bin.NtPerDay+avArgs.avRan*bin.NtPerDay-1;
%             bin.uPert(:,:,intervalBoundary)=0;
%             bin.vPert(:,:,intervalBoundary)=0;
%         end
                
    % Use RSCAT running mean method
    elseif avArgs.runMean==2
        
       % Check running mean not too long.
        if avArgs.avRan*2>=Nt
            error('Running mean longer than number of time cells');
        end
        
        if ~isfield(bin,'nonZeroIndex')
            fprintf('Calculating closest neighbours. \n');
            [bin.nonZeroIndex, bin.closestNeighbour]=neighbours(bin.count,...
                bin.NtPerDay,avArgs.avRan,kInd,bin.dateCell);
        end
        
        uRanAv=nan(size(bin.u));
        vRanAv=nan(size(bin.u));
        
        for m=1:size(bin.nonZeroIndex,1)
            i=bin.nonZeroIndex(m,1);
            j=bin.nonZeroIndex(m,2);
            k=bin.nonZeroIndex(m,3);
            uRanAv(i,j,k)=mean(bin.u(i,j,bin.closestNeighbour(m,:)));
            vRanAv(i,j,k)=mean(bin.v(i,j,bin.closestNeighbour(m,:)));
        end

        bin.uPert=bin.u-uRanAv;
        bin.vPert=bin.v-vRanAv;
       
        bin.uPert(isinf(bin.uPert) | isnan(bin.uPert) | bin.count==0)=0;
        bin.vPert(isinf(bin.vPert) | isnan(bin.vPert) | bin.count==0)=0;
              
    end
    
    clear uSum vSum uRanAv vRanAv countRanSum countRanSumAlt...
        countRanSumMax intervalBoundary
        
    fprintf('Calculating perturbation composites. \n');

    bin.uPertComp=zeros(length(bin.x),length(bin.y),bin.NtPerDay);
    bin.vPertComp=zeros(length(bin.x),length(bin.y),bin.NtPerDay);
    countPertCompAlt=zeros(length(bin.x),length(bin.y),bin.NtPerDay);

    for i=1:bin.NtPerDay
        bin.uPertComp(:,:,i)=sum(bin.uPert(:,:,i+(0:bin.numDaysTot-1)*bin.NtPerDay),3);
        bin.vPertComp(:,:,i)=sum(bin.vPert(:,:,i+(0:bin.numDaysTot-1)*bin.NtPerDay),3);
        % Count the places where uPert is non zero.
        countPertCompAlt(:,:,i)=sum(bin.uPert(:,:,i+(0:bin.numDaysTot-1)*bin.NtPerDay)~=0,3);
    end
    
    bin.uPertComp=bin.uPertComp./countPertCompAlt;
    bin.vPertComp=bin.vPertComp./countPertCompAlt;
    
    clear uSum vSum countPertComp countPert countPertCompAlt...
    countPertCompMax

    bin.uPertComp(isnan(bin.uPertComp) | isinf(bin.uPertComp))=0;
    bin.vPertComp(isnan(bin.vPertComp) | isinf(bin.uPertComp))=0;
    
    fprintf('Calculating perturbation variances. \n');

    uPertRes=bin.uPert-repmat(bin.uPertComp,[1 1 bin.numDaysTot]);
    vPertRes=bin.vPert-repmat(bin.vPertComp,[1 1 bin.numDaysTot]);

    % Remove residuals where we are simply subtracting binAv.uComp from zero
    uPertRes(bin.count==0)=0;
    vPertRes(bin.count==0)=0;
    
    % Calculate sample variance and covariance
    for i=1:bin.NtPerDay
        bin.uPertVar(:,:,i)=sum(uPertRes(:,:,i+(0:bin.numDaysTot-1)*bin.NtPerDay).^2,3);
        bin.vPertVar(:,:,i)=sum(vPertRes(:,:,i+(0:bin.numDaysTot-1)*bin.NtPerDay).^2,3);
        bin.pertCoVar(:,:,i)=sum(uPertRes(:,:,i+(0:bin.numDaysTot-1)*bin.NtPerDay).*...
            vPertRes(:,:,i+(0:bin.numDaysTot-1)*bin.NtPerDay),3);
    end

    bin.uPertVar=bin.uPertVar./(bin.countComp-1);
    bin.vPertVar=bin.vPertVar./(bin.countComp-1);
    bin.pertCoVar=bin.pertCoVar./(bin.countComp-1); 

    bin.uPertVar(isnan(bin.uPertVar) | isinf(bin.uPertVar) | bin.countComp<=1)=0;
    bin.vPertVar(isnan(bin.vPertVar) | isinf(bin.vPertVar) | bin.countComp<=1)=0;
    bin.pertCoVar(isnan(bin.pertCoVar) | isinf(bin.pertCoVar) | bin.countComp<=1)=0;
   
    clear uPertRes vPertRes
    
end

% Calculate T2 test statistic if required
if avArgs.t2Test
    
    fprintf('Calculating T2 test statistic. \n');
    
    if length(kInd)>2
        warning(['More than two time bins with data.'...
            ' Setting pValue to 0 \n']);
        bin.pValue=zeros(length(bin.x),length(bin.y));
    else
        bin.pValue=hotelling(bin,bin,avArgs.testTime1,avArgs.testTime2);
    end    
% Otherwise calculate manova test!
elseif avArgs.manTest    
    fprintf('Calculating Manova test statistic. \n');
    
    bin.pValue=nan(length(bin.x),length(bin.y));
    bin.pValue=multiMan(bin);
    
end

if isfield(bin,'pPert')
    bin=rmfield(bin,{'pPert'});
end
bin=rmfield(bin,{'u','v','count'});

if isfield(bin,'uPert')
    bin=rmfield(bin,{'uPert','vPert'});
end

% Calculate divergence if required.
if avArgs.divergence
    
    fprintf('Calculating divergence. \n');
    
    dLat=bin.y(2)-bin.y(1);
    dLon=bin.x(2)-bin.x(1);
    
    bin.uComp(bin.uComp==0)=NaN;
    bin.vComp(bin.vComp==0)=NaN;
   
    vPlus=cat(2,bin.vComp(:,2:length(bin.y),:),zeros(length(bin.x),1,bin.NtPerDay));
    vMinus=cat(2,zeros(length(bin.x),1,bin.NtPerDay),bin.vComp(:,1:length(bin.y)-1,:));
    uPlus=cat(1,bin.uComp(2:length(bin.x),:,:),zeros(1,length(bin.y),bin.NtPerDay));
    uMinus=cat(1,ones(1,length(bin.y),bin.NtPerDay),bin.uComp(1:length(bin.x)-1,:,:));
    
    bin.div=(vPlus-vMinus)./(2*dLat*110.57*10^3)+(uPlus-uMinus)./(2*dLon*111.32*10^3);
    
    clear uPlus uMinus vPlus vMinus
   
end

bin.avArgs=avArgs;

end
