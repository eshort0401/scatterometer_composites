function [mean, comp, pert, lat, lon] = binERA

% OLR appears to reset in 12 hour blocks with ERA-interim

% olrEnd stores the acculated ttr from the previous time step.
olrEnd=ncread('eraiOLRstatic.nc','ttr');
olrEnd=olrEnd(:,:,end);
olrEnd=repmat(olrEnd,[1 1 1]);
lat=ncread('eraiOLR.nc','latitude'); 
lon=ncread('eraiOLR.nc','longitude');

% olr stores the accumulated ttr 
olr=ncread('eraiOLR.nc','ttr');

% Load time array
t=ncread('eraiOLR.nc','time');

olrDaily=nan(size(olr,1),size(olr,2),size(olr,3)/8);
tDaily=nan(1,size(olr,3)/8);

for i=1:size(olr,3)/8
    olrDaily(:,:,i)=olr(:,:,(i-1)*8+4)+olr(:,:,(i-1)*8+8);
    tDaily(i)=nanmean(t((i-1)*8+(1:8)));
end 

olrDaily=olrDaily./(24*3600);

MJOdata=csvread('MJOdata.txt');

ascatStart=[1990 01 01 0 0 0];

tMJO=24*floor((tDaily-etime(ascatStart,[1900 01 01 0 0 0])/(60*60))/24);

mjoPhase=double(MJOdata((tMJO./24)+1,4));
mjoAmp=double(MJOdata((tMJO./24)+1,5));

tAus1=tMJO>=etime([2012 11 01 0 0 0],ascatStart)/(3600) & tMJO<etime([2013 05 01 0 0 0],ascatStart)/(3600);
tAus2=tMJO>=etime([2013 11 01 0 0 0],ascatStart)/(3600) & tMJO<etime([2014 05 01 0 0 0],ascatStart)/(3600);
tAus3=tMJO>=etime([2014 11 01 0 0 0],ascatStart)/(3600) & tMJO<etime([2015 05 01 0 0 0],ascatStart)/(3600);
tAus4=tMJO>=etime([2015 11 01 0 0 0],ascatStart)/(3600) & tMJO<etime([2016 05 01 0 0 0],ascatStart)/(3600);
tAus5=tMJO>=etime([2016 11 01 0 0 0],ascatStart)/(3600) & tMJO<etime([2017 05 01 0 0 0],ascatStart)/(3600);

tAus=tAus1 | tAus2 | tAus3 | tAus4 | tAus5;

tAs1=tMJO>=etime([2013 05 01 0 0 0],ascatStart)/(3600) & tMJO<etime([2013 11 01 0 0 0],ascatStart)/(3600);
tAs2=tMJO>=etime([2014 05 01 0 0 0],ascatStart)/(3600) & tMJO<etime([2014 11 01 0 0 0],ascatStart)/(3600);
tAs3=tMJO>=etime([2015 05 01 0 0 0],ascatStart)/(3600) & tMJO<etime([2015 11 01 0 0 0],ascatStart)/(3600);
tAs4=tMJO>=etime([2016 05 01 0 0 0],ascatStart)/(3600) & tMJO<etime([2016 11 01 0 0 0],ascatStart)/(3600);

tAs=tAs1 | tAs2 | tAs3 | tAs4;

sig=(mjoAmp>=1);

active1=olrDaily;
inactive1=olrDaily;

m1=(mjoPhase==3 | mjoPhase==4 | mjoPhase==5);

active1(lon>=94 & lon<110,:,~(m1 & sig))=NaN;
inactive1(lon>=94 & lon<110,:,~(~m1 & sig))=NaN;

active2=olrDaily;
inactive2=olrDaily;

m2=(mjoPhase==4 | mjoPhase==5);

active2(lon>=110 & lon<120,:,~(m2 & sig))=NaN;
inactive2(lon>=110 & lon<120,:,~(~m2 & sig))=NaN;

active3=olrDaily;
inactive3=olrDaily;

m3=(mjoPhase==4 | mjoPhase==5 | mjoPhase==6);

active3(lon>=120 & lon<140,:,~(m3 & sig))=NaN;
inactive3(lon>=120 & lon<140,:,~(~m3 & sig))=NaN;

active4=olrDaily;
inactive4=olrDaily;

m4=(mjoPhase==4 | mjoPhase==5 | mjoPhase==6 | mjoPhase==7);

active4(lon>=140,:,~(m4 & sig))=NaN;
inactive4(lon>=140,~(~m4 & sig))=NaN;

active=olrDaily;
inactive=olrDaily;

active(lon>=94 & lon<110,:,:)=active1(lon>=94 & lon<110,:,:);
active(lon>=110 & lon<120,:,:)=active2(lon>=110 & lon<120,:,:);
active(lon>=120 & lon<140,:,:)=active3(lon>=120 & lon<140,:,:);
active(lon>=140,:,:)=active4(lon>=140,:,:);

inactive(lon>=94 & lon<110,:,:)=inactive1(lon>=94 & lon<110,:,:);
inactive(lon>=110 & lon<120,:,:)=inactive2(lon>=110 & lon<120,:,:);
inactive(lon>=120 & lon<140,:,:)=inactive3(lon>=120 & lon<140,:,:);
inactive(lon>=140,:,:)=inactive4(lon>=140,:,:);

activeAus=active;
activeAs=active;
inactiveAus=inactive;
inactiveAs=inactive;

activeAus(:,:,~tAus)=NaN;
inactiveAus(:,:,~tAus)=NaN;
activeAs(:,:,~tAs)=NaN;
inactiveAs(:,:,~tAs)=NaN;

activeAus=nanmean(activeAus,3);
inactiveAus=nanmean(inactiveAus,3);
activeAs=nanmean(activeAs,3);
inactiveAs=nanmean(inactiveAs,3);

mean=nanmean(olrDaily,3);

meanAus=olrDaily;
meanAus(:,:,~tAus)=NaN;
meanAus=nanmean(meanAus,3);

meanAs=olrDaily;
meanAs(:,:,~tAs)=NaN;
meanAs=nanmean(meanAs,3);

activeAusPert=activeAus-mean;
inactiveAusPert=inactiveAus-mean;
activeAsPert=activeAs-mean;
inactiveAsPert=inactiveAs-mean;

mean={mean,meanAus,meanAs};
pert={activeAusPert,inactiveAusPert,activeAsPert,inactiveAsPert};
comp={activeAus,inactiveAus,activeAs,inactiveAs};

end

