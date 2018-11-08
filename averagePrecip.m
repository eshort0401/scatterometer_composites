function bin = averagePrecip(bin,avArgs)

%--------------------------------------------------------------------------
% Calculate composites and statistics
%--------------------------------------------------------------------------

if nargin<2
    avArgs.composite=1;
    avArgs.pertComp=1;
    
    % How many days either side to use for running mean.
    avArgs.avRan=.5;
    avArgs.runMean=0;
end

Nt=bin.numDaysTot*bin.NtPerDay;

fprintf('Calculating composites. \n');

% Calculate composites.
bin.pComp=zeros(length(bin.x),length(bin.y),bin.NtPerDay);
bin.tComp=zeros(length(bin.x),length(bin.y),bin.NtPerDay);

% bin.tVar=zeros(length(bin.x),length(bin.y),bin.NtPerDay);
bin.pVar=zeros(length(bin.x),length(bin.y),bin.NtPerDay);
bin.tVar=zeros(length(bin.x),length(bin.y),bin.NtPerDay);

% Create array to store composite count.    
bin.countComp=zeros(length(bin.x),length(bin.y),bin.NtPerDay);
if ~isfield(bin,'time')
    bin.time=zeros(length(bin.x),length(bin.y),Nt);
end

for i=1:bin.NtPerDay
    bin.pComp(:,:,i)=sum(bin.p(:,:,i+(0:bin.numDaysTot-1)*bin.NtPerDay),3);
    bin.tComp(:,:,i)=sum(bin.time(:,:,i+(0:bin.numDaysTot-1)*bin.NtPerDay),3);
    bin.countComp(:,:,i)=sum(bin.count(:,:,i+(0:bin.numDaysTot-1)*bin.NtPerDay)>0,3);
end

bin.pComp=bin.pComp./bin.countComp;
bin.tComp=bin.tComp./bin.countComp;

bin.pComp(isnan(bin.pComp) | isinf(bin.pComp))=0;
bin.tComp(isnan(bin.tComp) | isinf(bin.tComp))=0;

fprintf('Calculating variances. \n');

pRes=bin.p-repmat(bin.pComp,[1 1 bin.numDaysTot]);
tRes=bin.time-repmat(bin.tComp,[1 1 bin.numDaysTot]);

% Remove residuals where we are simply subtracting binAv.uComp from zero
pRes(bin.count==0)=0;
tRes(bin.count==0)=0;

% Calculate sample variance and covariance
for i=1:bin.NtPerDay
    bin.pVar(:,:,i)=sum(pRes(:,:,i+(0:bin.numDaysTot-1)*bin.NtPerDay).^2,3);
    bin.tVar(:,:,i)=sum(tRes(:,:,i+(0:bin.numDaysTot-1)*bin.NtPerDay).^2,3);
end

bin.pVar=bin.pVar./(bin.countComp-1);
bin.tVar=bin.tVar./(bin.countComp-1); 

bin.pVar(isnan(bin.pVar) | isinf(bin.pVar) | bin.countComp<=1)=0;
bin.tVar(isnan(bin.tVar) | isinf(bin.tVar) | bin.countComp<=1)=0;

clear pRes tRes

bin=rmfield(bin,'time');
startTime=bin.dateCell{1};
bin.tComp(bin.tComp==0)=NaN;
bin.tComp=mod(bin.tComp+startTime(4)+startTime(5)/60+startTime(6)/3600,24);

%--------------------------------------------------------------------------
% Calculate perturbations if required
%--------------------------------------------------------------------------

% Calculate perturbation winds if required. Be careful to not double weight
% the current time. For now ignore future most time. 
if avArgs.pertComp
    
    fprintf('Calculating perturbation winds. \n');
    
    kInd=find(squeeze(squeeze(any(any(bin.countComp>0,1),2)))==1);
    kInd=unique(kInd);
   
    bin.pPert=zeros(length(bin.x),length(bin.y),Nt);
                
    % Use basic mean.  
    if avArgs.runMean==0
        
        pCompNan=bin.pComp;
        
        % Check there is data at each time.
        pCompNan(bin.pComp==0)=NaN;

        pRanAv=nanmean(pCompNan,3);

        clear uCompNan vCompNan
        
        for i=1:size(bin.p,3)
            bin.pPert(:,:,i)=bin.p(:,:,i)-pRanAv;
        end
        
        bin.pPert(isnan(bin.pPert) | isinf(bin.pPert) | bin.count==0)=0;
        bin.vPert(isnan(bin.pPert) | isinf(bin.pPert) | bin.count==0)=0;
       
    % Use simple running mean method.
    elseif avArgs.runMean==1
        % Check running mean not too long.
        if avArgs.avRan*2>=Nt
            error('Running mean longer than number of time cells');
        end

        for k=avArgs.avRan*bin.NtPerDay+1:Nt-avArgs.avRan*bin.NtPerDay+1
            
            countRanAlt=sum(bin.count(:,:,k-avArgs.avRan*bin.NtPerDay:k+avArgs.avRan*bin.NtPerDay-1)>0,3);

            pRanAv=sum(bin.u(:,:,k-avArgs.avRan*bin.NtPerDay:k+avArgs.avRan*bin.NtPerDay-1),3);
            pRanAv=pRanAv./countRanAlt;
            
            vRanAv=sum(bin.v(:,:,k-avArgs.avRan*bin.NtPerDay:k+avArgs.avRan*bin.NtPerDay-1),3);
            vRanAv=vRanAv./countRanAlt;
            
            uPertK=bin.u(:,:,k)-pRanAv;
            vPertK=bin.v(:,:,k)-vRanAv;
            
            uPertK(countRanAlt<avArgs.runThreshAlt)=0;
            vPertK(countRanAlt<avArgs.runThreshAlt)=0;
            
            uPert(:,:,k)=uPertK;
            vPert(:,:,k)=vPertK;
            
        end
        
        uPert(isinf(uPert) | isnan(uPert) | bin.count==0)=0;
        vPert(isinf(uPert) | isnan(uPert) | bin.count==0)=0;
        
        % Remove runMean values around the ends of disjoint intervals 
        for i=1:length(bin.numDays)-1
            intervalBoundary=sum(bin.numDays(1:i))*bin.NtPerDay-...
                avArgs.avRan*bin.NtPerDay:sum(bin.numDays(1:i))*...
                bin.NtPerDay+avArgs.avRan*bin.NtPerDay-1;
            uPert(:,:,intervalBoundary)=0;
            vPert(:,:,intervalBoundary)=0;
        end
                
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
        
        pRanAv=nan(size(bin.p));
        
        for m=1:size(bin.nonZeroIndex,1)
            i=bin.nonZeroIndex(m,1);
            j=bin.nonZeroIndex(m,2);
            k=bin.nonZeroIndex(m,3);
            pRanAv(i,j,k)=mean(bin.p(i,j,bin.closestNeighbour(m,:)));
        end

        bin.pPert=bin.p-pRanAv;
        
        bin.pPert(isinf(bin.pPert) | isnan(bin.pPert) | bin.count==0)=0;
              
    end
    
    clear pRanAv countRanSum countRanSumAlt...
        countRanSumMax intervalBoundary
        
    fprintf('Calculating perturbation composites. \n');

    bin.pPertComp=zeros(length(bin.x),length(bin.y),bin.NtPerDay);
    countPertCompAlt=zeros(length(bin.x),length(bin.y),bin.NtPerDay);

    for i=1:bin.NtPerDay
        bin.pPertComp(:,:,i)=sum(bin.pPert(:,:,i+(0:bin.numDaysTot-1)*bin.NtPerDay),3);
        % Count the places where uPert is non zero.
        countPertCompAlt(:,:,i)=sum(bin.count(:,:,i+(0:bin.numDaysTot-1)*bin.NtPerDay)>0,3);
    end
    
    bin.pPertComp=bin.pPertComp./countPertCompAlt;
    
    clear uSum vSum countPertComp countPert countPertCompAlt...
    countPertCompMax

    bin.pPertComp(isnan(bin.pPertComp) | isinf(bin.pPertComp))=0;
    
    fprintf('Calculating perturbation variances. \n');

    pPertRes=bin.pPert-repmat(bin.pPertComp,[1 1 bin.numDaysTot]);

    % Remove residuals where we are simply subtracting binAv.uComp from zero
    pPertRes(bin.count==0)=0;
    
    % Calculate sample variance and covariance
    for i=1:bin.NtPerDay
        bin.pPertVar(:,:,i)=sum(pPertRes(:,:,i+(0:bin.numDaysTot-1)*bin.NtPerDay).^2,3);
    end

    bin.pPertVar=bin.pPertVar./(bin.countComp-1);

    bin.pPertVar(isnan(bin.pPertVar) | isinf(bin.pPertVar) | bin.countComp<=1)=0;
   
    clear pPertRes
    
end

bin=rmfield(bin,{'p','count'});
bin.pValue=zeros(length(bin.x),length(bin.y));
bin.avArgs=avArgs;

end

